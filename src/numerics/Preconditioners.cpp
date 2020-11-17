#include "cantera/numerics/Preconditioners.h"

namespace Cantera
{
    /**
     * 
     * Preconditioner Implementation
     * 
     * **/
    double PreconditionerBase::getThreshold()
    {
        return this->threshold;
    }

    void PreconditionerBase::setThreshold(double threshold)
    {
        this->threshold = threshold;
    }

    void PreconditionerBase::setElementByThreshold(unsigned long row,unsigned long col, double element)
    {
        if (this->threshold<element)
        {
            this->setElement(row,col,element);
        }
    }

    void PreconditionerBase::setDimensions(unsigned long nrows,unsigned long ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
    }

    unsigned long* PreconditionerBase::getDimensions()
    {
        return &(this->dimensions[0]);
    }
}


namespace Cantera::AMP //Making ASP apart of Cantera namespace
{
    
    /**
     * 
     * AdaptivePreconditioner implementations
     * 
     * **/

    AdaptivePreconditioner::AdaptivePreconditioner()
    {
        this->addToFunctionMap("temperature",TemperatureDerivatives); //Adding temperature to function map
    }

    void AdaptivePreconditioner::setDimensions(unsigned long nrows,unsigned long ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
    }

    void AdaptivePreconditioner::setElement(unsigned long row,unsigned long col,double element)
    {
        this->matrix.coeffRef(row,col)=element;
    }

    double AdaptivePreconditioner::getElement(unsigned long row,unsigned long col)
    {
        warn_user("AdaptivePreconditioner::getElement","getElement not properly implemented yet, returning 1.0");
        return 1.0;//FIXME
    }

    Eigen::SparseMatrix<double>* AdaptivePreconditioner::getMatrix()
    {
        return &(this->matrix);
    }

    void AdaptivePreconditioner::setMatrix(Eigen::SparseMatrix<double> *sparseMat)
    {
        this->matrix=*(sparseMat);
    }

    void AdaptivePreconditioner::solve(double* x, double* b,unsigned long size)
    {   
        //b is currently the mass fractions that come in 
        
        //Compressing sparse matrix structure
        this->matrix.makeCompressed();
        //Creating vectors in the form of //Ax=b
        Eigen::Map<Eigen::VectorXd> bVector(b,size);
        // Eigen::Map<Eigen::VectorXd> xVector(x,size);
        Eigen::VectorXd xVector;
        //Create eigen solver object
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        // Initialize solver object
        solver.compute(this->matrix);
        Cantera::AMP::checkEigenError("AdaptivePreconditioner::solve (Eigen Decomposition)",solver.info());
        //Solve for xVector
        xVector = solver.solve(bVector);
        Cantera::AMP::checkEigenError("AdaptivePreconditioner::solve (Eigen Solve)",solver.info());
        //Copy x vector to x
        Eigen::VectorXd::Map(x,xVector.rows())= xVector;
    }

    void AdaptivePreconditioner::setup(Reactor *reactor, double t, double* y, double* ydot, double* params, unsigned long reactorStart)
    {   
        double *inputs[4] = {&t,y,ydot,params}; //Array of input arrays
        StateMap stateMap = getStateMap(reactor,reactorStart);
        // printReactorComponents(reactor);
       //Getting species on species derivatives
       Cantera::AMP::SpeciesDerivatives(this,reactor,inputs,stateMap);    
        //Solving other variables
        for (unsigned long i = 0; i < stateMap["species"]; i++)
        {
            //Getting component name
            std::string component = reactor->componentName(i);
            //If key not found in function map, no precondition
            if ( this->functionMap.find(component) == this->functionMap.end() ) 
            {
                // std::cout<<"No Precondition: "<<component<<std::endl;
                NoPrecondition(this,reactor,inputs,stateMap,component);
            } 
            //Key is found call appropriate preconditioner function
            else {
                //Calling component function
                this->functionMap[component](this,reactor,inputs,stateMap,component);
            }  
        }
        
        std::cout<<Eigen::MatrixXd(this->matrix)<<std::endl;   
    }

    void AdaptivePreconditioner::initialize(unsigned long nrows,unsigned long ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
        this->nonzeros = this->dimensions[0]*this->dimensions[1]/2; //Reserves up to half the total spaces
        this->matrix.resize(nrows,ncols);
        this->matrix.reserve(this->nonzeros);
    }

    void AdaptivePreconditioner::reset()
    {
        //Do any reset stuff here
        this->matrix.setZero(); //Set all elements to zero
        this->matrix.makeCompressed(); //Compress matrix
        this->matrix.reserve(this->nonzeros); //Reserve space potentially needed
    }

    void AdaptivePreconditioner::addToFunctionMap(std::string component, AdaptiveFunction newFunction)
    {
        this->functionMap[component] = newFunction;
    }

    void AdaptivePreconditioner::removeFromFunctionMap(std::string component)
    {
        if ( this->functionMap.find(component) != this->functionMap.end() ) 
        {
            this->functionMap.erase(component);
        }
        else
        {
            warn_user("Cantera::AdaptivePreconditioner::removeFromFunctionMap",component+": key not found");
        }
    }

    /**
     * 
     * Non-Class member functions
     * 
     * */
    void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap stateMap, std::string key)
    {   
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        ThermoPhase* thermo=reactor->getThermoMgr();
        //Important sizes to the determination of values
        unsigned long numberOfSpecies = kinetics->nTotalSpecies();
        unsigned long speciesStart = stateMap["species"];
        unsigned long tempIndex = stateMap["temperature"]-1;

        //Array pointers for data that is reused
        //net production rates (omega dot)
        double* netProductionRatesNext = new double[numberOfSpecies];
        double* netProductionRatesCurrent = new double[numberOfSpecies];
        
        //Getting perturbed state
        //Perturbation for finite difference of temperature
        double deltaTemp = inputs[1][tempIndex]*(std::sqrt(__DBL_EPSILON__));
        thermo->setTemperature(inputs[1][tempIndex]+deltaTemp);
        kinetics->getNetProductionRates(netProductionRatesNext);
        double TDotNext = reactor->evaluateEnergyEquation(*inputs[0],inputs[1],inputs[2],inputs[3])/(thermo->cp_mass()*reactor->mass()); //Perturbed internal energy

        //Getting current state
        thermo->setTemperature(inputs[1][tempIndex]); //Setting temperature back to correct value
        kinetics->getNetProductionRates(netProductionRatesCurrent);
        double TDotCurrent = reactor->evaluateEnergyEquation(*inputs[0],inputs[1],inputs[2],inputs[3])/(thermo->cp_mass()*reactor->mass()); //Current internal energy

        /**
         * Temp Rate Derivatives w.r.t Temp
         * d T_dot/dT
         **/
        temperatureTemperatureDerivatives(preconditioner,TDotCurrent,TDotNext,deltaTemp,tempIndex);

        /**
         * Production Rate Derivatives w.r.t Temp
         * d omega_dot_j/dT
         **/
        //convert kmol/m^3/s to kmol/s by multiplying volume and deltaTemp
        speciesTemperatureDerivatives(preconditioner,netProductionRatesCurrent,netProductionRatesNext,deltaTemp*reactor->volume(),numberOfSpecies,speciesStart,tempIndex);
        
        //Deleting appropriate pointers
        delete[] netProductionRatesNext;
        delete[] netProductionRatesCurrent;
    }
        
    inline void temperatureTemperatureDerivatives(PreconditionerBase *preconditioner, double TDotCurrent, double TDotNext, double deltaTemp, unsigned long tempIndex)
    {
        preconditioner->setElementByThreshold(tempIndex,tempIndex,(TDotNext-TDotCurrent)/deltaTemp); 
    }


    inline void speciesTemperatureDerivatives(PreconditionerBase *preconditioner, double *netProductionRatesCurrent, double *netProductionRatesNext, double deltaTemp,unsigned long numberOfSpecies,unsigned long speciesStart, unsigned long tempIndex)
    {
        for (unsigned long j = 0; j < numberOfSpecies; j++) //column
        {   
            preconditioner->setElementByThreshold(j+speciesStart,tempIndex,(netProductionRatesNext[j]-netProductionRatesCurrent[j])/deltaTemp); //Add by threshold specTempDerivative
        }

    }


    /*
        Species Derivatives
    */
    void SpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs, StateMap stateMap)
    {
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        //Getting thermophase object for access to concentrations and species data
        ThermoPhase* thermo=reactor->getThermoMgr();
        //Compositions for reactants and products
        Composition reactants;
        Composition products;
        //Important sizes to the determination of values
        unsigned long numberOfReactions = kinetics->nReactions();
        unsigned long numberOfSpecies = kinetics->nTotalSpecies();     
        //Array pointers for data that is reused
        double* rateLawDerivatives = new double[numberOfSpecies*numberOfSpecies];
        double* kForward = new double[numberOfReactions];
        double* kBackward = new double[numberOfReactions];
        //Concentrations of species
        double* concentrations = new double[numberOfSpecies];
        //Getting species concentrations
        thermo->getConcentrations(concentrations);
        //Getting forward rate constants for calcs
        kinetics->getFwdRateConstants(kForward); 
        //Getting reverse rate constants for calcs
        kinetics->getRevRateConstants(kBackward); 
        //shared_ptr for current reaction in finding Jacobian
        std::shared_ptr<Reaction> currentReaction;
        for (unsigned long r = 0; r < numberOfReactions; r++)
        {
            currentReaction=kinetics->getReactionPtr(r);
            //Loop through reactants in current reaction
            reactants = currentReaction->reactants;
            products = currentReaction->products;
            Cantera::AMP::speciesSpeciesDerivative(reactants,stateMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume(),numberOfSpecies);
            Cantera::AMP::speciesSpeciesDerivative(products,stateMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume(),numberOfSpecies); //Multiply by negative one to change direction
        }

        temperatureSpeciesDerivatives(preconditioner,concentrations,rateLawDerivatives,reactor->volume(),thermo,kinetics,stateMap);
        //Adding to preconditioner
        //d(w)/dn_j
        unsigned long idx;
        unsigned long speciesStart = stateMap["start"]+stateMap["species"];
        for (unsigned long j = 0; j < numberOfSpecies; j++) // column
            {
            for (unsigned long i = 0; i < numberOfSpecies; i++) //row
            {  
            idx = j+i*numberOfSpecies; //Getting flattened index
            preconditioner->setElementByThreshold(i+speciesStart,j+speciesStart,reactor->volume()*rateLawDerivatives[idx]);//Add by threshold 
            }
        }

        //Deleting appropriate pointers
        delete[] kForward;
        delete[] kBackward;
        delete[] concentrations;
        delete[] rateLawDerivatives;
    }

    inline void temperatureSpeciesDerivatives(PreconditionerBase *preconditioner, double *concentrations, double *dwkdnj, double volume, ThermoPhase *thermo,Kinetics *kinetics, StateMap stateMap)
    {
        /**
         * Temp Rate Derivatives w.r.t Species
         * d T_dot/dnj
         **/
        unsigned long numberOfSpecies = thermo->nSpecies();
        unsigned long tempIndex = stateMap["temperature"]-1;
        unsigned long speciesStart = stateMap["species"];
        double *enthalpy = new double[numberOfSpecies];
        double *specificHeat = new double[numberOfSpecies];
        double *netProductionRates = new double[numberOfSpecies];
        thermo->getPartialMolarEnthalpies(enthalpy);
        thermo->getPartialMolarCp(specificHeat);
        kinetics->getNetProductionRates(netProductionRates);
        //Getting perturbed changes w.r.t temperature
        for (unsigned long j = 0; j < numberOfSpecies; j++) //Spans columns
        {
            double CkCpkSum = 0;
            double hkwkSum = 0;
            double hkdwkdnjSum = 0;
            for (unsigned long k = 0; k < numberOfSpecies; k++) //Spans rows
            {   
                int idx = j+k*numberOfSpecies; //Getting flattened index - j to remain same and k to change. This means moving down a row.
                hkwkSum += enthalpy[k]*netProductionRates[k];
                hkdwkdnjSum += enthalpy[k]*dwkdnj[idx];
                CkCpkSum += concentrations[k]*specificHeat[k];
            }
            //Set appropriate colume of preconditioner
            preconditioner->setElementByThreshold(tempIndex,j+speciesStart,(-hkdwkdnjSum*CkCpkSum+specificHeat[j]/volume*hkwkSum)/(CkCpkSum*CkCpkSum));
        }
        delete[] enthalpy;
        delete[] specificHeat;
        delete[] netProductionRates;
    }

    inline void speciesSpeciesDerivative(std::map<std::string, double> comp,StateMap stateMap, double* omega, double* concentrations, double k_direction, double volume, unsigned long numberOfSpecies)
    { 
        //flattened index for derivatives
        unsigned long oidx; //index for omega
        unsigned long sidx; //index for species
        unsigned long speciesStart = stateMap["species"]; //Adjustment based on reactor, network, and stateMap to get Omega index
        double dRdn; //temporary value for rate derivative
        for (std::map<std::string,double>::iterator iter1 = comp.begin(); iter1 != comp.end(); iter1++) //Independent variable -- column
        {
            for (std::map<std::string,double>::iterator iter2 = comp.begin(); iter2 != comp.end(); iter2++) //Dependent variable -- row
            {
            //Get index for current omega
            oidx = stateMap[iter1->first]-speciesStart+(stateMap[iter2->first]-speciesStart)*numberOfSpecies;
            dRdn=1; //Set dRdn to one so multiplication isn't zero unless a coefficient is zero
            // Loop through species to get rate derivative
            for (std::map<std::string,double>::iterator iter3 = comp.begin(); iter3 != comp.end(); iter3++) //Dependent variable
            {
                sidx = stateMap[iter3->first]-speciesStart;
                if (iter3->first == iter1->first)
                {
                dRdn *= std::pow(iter3->second*concentrations[sidx],iter3->second-1); //derivative
                }
                else
                {
                dRdn *= std::pow(concentrations[sidx],iter3->second); //Not derivative
                }
            }
            omega[oidx] += k_direction*dRdn/volume; //Updating omega derivative as is necessary
            }
        }
    } 

    

    void NoPrecondition(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs, StateMap stateMap, std::string key)
    {   
        unsigned long idx = reactor->componentIndex(key)+stateMap["start"];
        preconditioner->setElement(idx,idx,1); //setting mass variable element of preconditioner equal to 1
    }

    /*
        Other functions used in preconditioner functions but not directly related to a state variable
    */

   void ForwardConversion(Reactor *currReactor, double *tempState, double *rhs, unsigned long m_start)
    {
        ThermoPhase *thermo = currReactor->getThermoMgr();
        double currMass =  currReactor->mass();
        unsigned long nStateVars = currReactor->neq()-thermo->nSpecies(); 
        //Transferring unchanged parameters for each reactor
        for (unsigned long i = 0; i < nStateVars; i++)
        {   
            unsigned long globalIndex = m_start+i;
            tempState[globalIndex] = rhs[globalIndex];
        }
        //Adjusting mass fraction parameters for each reactor to moles for AJP
        double *molecularWeights = new double[thermo->nSpecies()];
        thermo->getMolecularWeights(molecularWeights);
        for (unsigned long i = 0; i < thermo->nSpecies(); i++)
        {
            unsigned long globalIndex = m_start+i+nStateVars;
            tempState[globalIndex] = currMass*molecularWeights[i]*rhs[globalIndex];     
        }
        delete[] molecularWeights;
    }

    void BackwardConversion(Reactor *currReactor, double *output, unsigned long m_start)
    {
            ThermoPhase *thermo = currReactor->getThermoMgr();
            double currMass =  currReactor->mass();
            unsigned long nStateVars = currReactor->neq()-thermo->nSpecies(); 
            //Do nothing to unchanged parameters
            //Convert moles back to mass fractions
            double *molecularWeights = new double[thermo->nSpecies()];
            thermo->getMolecularWeights(molecularWeights);
            for (unsigned long i = 0; i < thermo->nSpecies(); i++)
            {
                unsigned long globalIndex = m_start+i+nStateVars;
                output[globalIndex] *= 1/(molecularWeights[i]*currMass);     
            }
            delete[] molecularWeights;
    }
    
    StateMap getStateMap(Reactor *reactor,unsigned long start)
    {
        StateMap stateMap;
        unsigned long nonspeciesNumber = reactor->neq()-reactor->getKineticsMgr()->nTotalSpecies();
        for (unsigned long i = 0; i < reactor->neq(); i++)
        {
            stateMap[reactor->componentName(i)] = i+1;
        }
        stateMap["species"]=nonspeciesNumber; //Start of species
        stateMap["start"] = start;
        return stateMap;
    }
    
    inline void printReactorComponents(Reactor* reactor)
    {
        for (unsigned long i = 0; i < reactor->neq(); i++)
        {
        std::cout<<reactor->componentName(i)<<std::endl;
        }
    }

    void checkEigenError(std::string method, unsigned long info)
    {
        if(info!=Eigen::Success) 
        {   
            std::string error="Failure: ";
            if(info==Eigen::NumericalIssue)
            {
                error+="NumericalIssues";
            }
            else if(info==Eigen::NoConvergence)
            {
                error+="NoConvergence";
            }
            else if(info==Eigen::InvalidInput)
            {
                error+="InvalidInput";
            }
            else
            {
                error+="Unknown";
            }
            warn_user(method,error);
            throw CanteraError(method,error);
        }
    }

}
