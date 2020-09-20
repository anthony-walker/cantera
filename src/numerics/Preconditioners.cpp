#include "cantera/numerics/Preconditioners.h"
#include <regex>

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

    void PreconditionerBase::setElement(unsigned long row, unsigned long col, double element)
    {
        warn_user("PreconditionerBase::setElement","setElement function has not been overriden");
    }


    unsigned long* PreconditionerBase::getDimensions()
    {
        return &(this->dimensions[0]);
    }

    double PreconditionerBase::getElement(unsigned long row, unsigned long col)
    {
        warn_user("PreconditionerBase::getElement","getElement function has not been overriden, returning 1.");
        return 1.0;
    }

    void PreconditionerBase::initialize(unsigned long nrows,unsigned long ncols)
    {
        warn_user("PreconditionerBase::initialize","initialize function has not been overriden.");
    }
    
    void PreconditionerBase::reset()
    {
        warn_user("PreconditionerBase::reset","reset function has not been overriden.");
    }
    
    void PreconditionerBase::setup(Reactor *reactor, double t, double* y, double* ydot, double* params, unsigned long reactorStart)
    {
        throw CanteraError("PreconditionerBase::setup", "Reactor type:(Reactor) is not implemented for the specified preconditioner type.");
    }

    void PreconditionerBase::solve(double* x, double* b,unsigned long size)
    {
        throw CanteraError("PreconditionerBase::setup", "Reactor type:(IdealGasConstPressureReactor) is not implemented for the specified preconditioner type.");
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
        this->functionMap["temperature"] = TemperatureDerivatives;
        this->functionMap["volume"]=NoPrecondition;
        this->functionMap["pressure"]=NoPrecondition;
        this->functionMap["mass"]=NoPrecondition;
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
        //rateLawDerivatives used in other functions
        unsigned long numberOfSpecies = reactor->getKineticsMgr()->nTotalSpecies();
        double* rateLawDerivatives = new double[numberOfSpecies*numberOfSpecies];
        IndexMap idxMap = getNonSpeciesIndexMap(reactor,reactorStart);
       //Getting species on species derivatives
       Cantera::AMP::SpeciesSpeciesDerivatives(this,reactor,y,ydot,rateLawDerivatives,idxMap);    
        //Solving other variables
        for (unsigned long i = 0; i < idxMap["species"]; i++)
        {
            //Getting component name
            std::string component = reactor->componentName(i);
            //Calling component function
            this->functionMap[component](this,reactor,y,ydot,rateLawDerivatives,idxMap,component);
        }
        //Deleting rateLawDerivatives array
        delete[] rateLawDerivatives;
        // std::cout<<Eigen::MatrixXd(this->matrix)<<std::endl;   
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

    /**
     * 
     * Non-Class member functions
     * 
     * */

    
    void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double* y, double* ydot, double* rateLawDerivatives,IndexMap indexMap, std::string key)
    {   
        // //Getting kinetics object for access to reactions
        // Kinetics* kinetics=reactor->getKineticsMgr();
        // ThermoPhase* thermo=reactor->getThermoMgr();
        // //Important sizes to the determination of values
        // unsigned long numberOfSpecies = kinetics->nTotalSpecies();

        // //Array pointers for data that is reused
        // //net production rates (omega dot)
        // double* netProductionRatesNext = new double[numberOfSpecies];
        // double* netProductionRatesCurrent = new double[numberOfSpecies];
        // double* molecularWeights = new double[numberOfSpecies];
        // double* ydotPerturbed = new double[reactor->neq()];
        // //Perturbation for finite difference of temperature
        // double deltaTemp = y[index]*(std::sqrt(__DBL_EPSILON__));
        // thermo->getMolecularWeights(molecularWeights);
        // //Getting current state
        // kinetics->getNetProductionRates(netProductionRatesCurrent);
        // double intEnergyCurrent = thermo->intEnergy_mass(); //Current internal energy
        // //Getting perturbed state
        // thermo->setTemperature(y[index]+deltaTemp);
        // kinetics->getNetProductionRates(netProductionRatesNext);
        // double intEnergyNext = thermo->intEnergy_mass(); //Perturbed internal energy
        // double inverseDensity = 1/thermo->density();
        // printf("%0.15f, %0.15f\n",intEnergyCurrent,intEnergyNext);
        // double energyTotal=0.0;
        // //Getting perturbed changes w.r.t temperature
        // for (unsigned long i = 0; i < numberOfSpecies; i++)
        // {
        //     ydotPerturbed[i+speciesStart]=netProductionRatesNext[i]*molecularWeights[i]*inverseDensity;
        //     energyTotal+=intEnergyNext/y[speciesStart+i]*netProductionRatesNext[i];
        // }
        // ydotPerturbed[index]=-Cantera::GasConstant*(y[index]+deltaTemp)*energyTotal*inverseDensity/thermo->cv_mass();
        // //Setting mass and volume of perturbed state to be same as current state so that they are unaffected
        // ydotPerturbed[0] = ydot[index-2];
        // ydotPerturbed[1] = ydot[index-1];


        // //Adding to preconditioner the species derivatives
        // //d(dTdt)/dn_j
        // for (unsigned long j = 0; j < numberOfSpecies; j++) //column
        // {   
        //     double tempSpecDerivative = 3; //+1 because specie index will start after temperature
        //     preconditioner->setElementByThreshold(index,j+speciesStart,tempSpecDerivative); //Add by threshold tempSpecDerivative
        // }
        // //d(n_j)/dT -- seems correct
        // for (unsigned long j = 0; j < numberOfSpecies; j++) //column
        // {   
        //     double specTempDerivative = (netProductionRatesNext[j]-netProductionRatesCurrent[j])/deltaTemp;
        //     preconditioner->setElementByThreshold(j+speciesStart,index,specTempDerivative); //Add by threshold specTempDerivative
        // }
        // //Adding to preconditioner the temperature derivative
        // dTdt -= ydot[index];
        // // dTdt /= perturbation; //Turning dTdt into d(dTdt)/dT
        // //d(dTdt)/dT
        // preconditioner->setElementByThreshold(index,index,1); //Add by threshold dTdt
        // //Setting temperature back to correct value
        // thermo->setTemperature(y[index]);
        // //Deleting appropriate pointers
        // delete[] netProductionRatesNext;
        // delete[] netProductionRatesCurrent;
        // delete[] ydotPerturbed;
        }

        void SpeciesSpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double* y, double* ydot, double* rateLawDerivatives,IndexMap indexMap)
        {
        // //Getting kinetics object for access to reactions
        // Kinetics* kinetics=reactor->getKineticsMgr();
        // //Getting thermophase object for access to concentrations and species data
        // ThermoPhase* thermo=reactor->getThermoMgr();
        // //Compositions for reactants and products
        // Composition reactants;
        // Composition products;
        // //Important sizes to the determination of values
        // unsigned long numberOfReactions = kinetics->nReactions();
        // unsigned long numberOfSpecies = kinetics->nTotalSpecies();
        // //Array pointers for data that is reused
        // double* kForward = new double[numberOfReactions];
        // double* kBackward = new double[numberOfReactions];
        // //Concentrations of species
        // double* concentrations = new double[numberOfSpecies];
        // //Getting species names
        // const std::vector<std::string> *names = &(thermo->speciesNames());
        // //Getting species concentrations
        // thermo->getConcentrations(concentrations);
        // //Getting forward rate constants for calcs
        // kinetics->getFwdRateConstants(kForward); 
        // //Getting reverse rate constants for calcs
        // kinetics->getRevRateConstants(kBackward); 
        // //shared_ptr for current reaction in finding Jacobian
        // std::shared_ptr<Reaction> currentReaction;
        // //Creating map of species indices
        // IndexMap indexMap;
        // for (unsigned long i = 0; i < numberOfSpecies; i++)
        // {
        //     indexMap[names->at(i)]=i;
        // }
        // for (unsigned long r = 0; r < numberOfReactions; r++)
        // {
        //     currentReaction=kinetics->getReactionPtr(r);
        //     //Loop through reactants in current reaction
        //     reactants = currentReaction->reactants;
        //     products = currentReaction->products;
        //     Cantera::AMP::speciesDerivative(reactants,indexMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume());
        //     Cantera::AMP::speciesDerivative(products,indexMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume()); //Multiply by negative one to change direction
        // }

        // //Adding to preconditioner
        // //d(w)/dn_j
        // unsigned long idx;
        // for (unsigned long j = 0; j < numberOfSpecies; j++) // column
        //     {
        //     for (unsigned long i = 0; i < numberOfSpecies; i++) //row
        //     {  
        //     idx = j+i*numberOfSpecies; //Getting flattened index
        //     preconditioner->setElementByThreshold(i+speciesStart,j+speciesStart,reactor->volume()*rateLawDerivatives[idx]);//Add by threshold 
        //     }
        // }

        // //Deleting appropriate pointers
        // delete[] kForward;
        // delete[] kBackward;
        // delete[] concentrations;
    }

    IndexMap getNonSpeciesIndexMap(Reactor *reactor,unsigned long start)
    {
        IndexMap idxMap;
        unsigned long nonspeciesNumber = reactor->neq()-reactor->getKineticsMgr()->nTotalSpecies();
        for (unsigned long i = 0; i < nonspeciesNumber; i++)
        {
            idxMap[reactor->componentName(i)] = i;
        }
        idxMap["species"]=nonspeciesNumber; //Start of species
        idxMap["start"] = start;
        return idxMap;
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

    inline void speciesDerivative(std::map<std::string, double> comp,IndexMap indexMap, double* omega, double* concentrations, double k_direction, double volume)
    { 
        //flattened index for derivatives
        unsigned long oidx; //index for omega
        unsigned long sidx; //index for species
        double dRdn; //temporary value for rate derivative

        for (std::map<std::string,double>::iterator iter1 = comp.begin(); iter1 != comp.end(); iter1++) //Independent variable -- column
        {
            for (std::map<std::string,double>::iterator iter2 = comp.begin(); iter2 != comp.end(); iter2++) //Dependent variable -- row
            {
            //Get index for current omega
            oidx = indexMap[iter1->first]+indexMap[iter2->first]*indexMap.size();
            dRdn=1; //Set dRdn to one so multiplication isn't zero unless a coefficient is zero
            //Loop through species to get rate derivative
            for (std::map<std::string,double>::iterator iter3 = comp.begin(); iter3 != comp.end(); iter3++) //Dependent variable
            {
                sidx = indexMap[iter3->first];
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

    void NoPrecondition(PreconditionerBase *preconditioner,Reactor* reactor, double* y, double* ydot, double* rateLawDerivatives,IndexMap indexMap, std::string key)
    {
        unsigned long idx = indexMap[key]+indexMap["start"];
        preconditioner->setElement(idx,idx,1); //setting mass variable element of preconditioner equal to 1
    }
}
