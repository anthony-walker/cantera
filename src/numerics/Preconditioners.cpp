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

    void PreconditionerBase::solve(double* x, double* b)
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

    void AdaptivePreconditioner::solve(double* x, double* b)
    {
        
        std::cout<<"Here1"<<std::endl;
        this->matrix.makeCompressed();
        //Creating vectors in the form of //Ax=b
        std::cout<<"Here2"<<std::endl;
        typedef Eigen::Vector<double,sizeof(b)/sizeof(b[0])> systemVector;
        systemVector bVector = Eigen::Map<systemVector>(b); //copy b to bVector
        Eigen::VectorXd xVector;
        //Create eigen solver object
        std::cout<<"Here4"<<std::endl;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        // Initialize solver object
        std::cout<<"Here4"<<std::endl;
        solver.compute(this->matrix);
        if(solver.info()!=Eigen::Success) 
        {
            throw CanteraError("AdaptivePreconditioner::solve","Eigen decomposition was unsuccessful");
        }
        //Solve for xVector
        std::cout<<"Here5"<<std::endl;
        std::cout << bVector[0] << std::endl;
        std::cout << xVector[0] << std::endl;
        xVector = solver.solve(bVector);
        if(solver.info()!=Eigen::Success) 
        {
            throw CanteraError("AdaptivePreconditioner::solve","Eigen solving linear system was unsuccessful");
        }
        //Solve fo
        // xVector = solver.solve(bVector);
        // //Copy x vector to x
        std::cout<<"Here6"<<std::endl;
        // Eigen::VectorXd::Map(x,xVector.rows()) = xVector;
        // std::cout<<"Here7"<<std::endl;
        // //reset preconditioner
        // this->reset(); 
        // std::cout<<"Here8"<<std::endl;
    }

    void AdaptivePreconditioner::setup(Reactor *reactor, double t, double* y, double* ydot, double* params, unsigned long reactorStart)
    {   
        unsigned long speciesStart = reactor->neq()-reactor->getKineticsMgr()->nTotalSpecies();
       //Getting species on species derivatives
       Cantera::AMP::SpeciesSpeciesDerivatives(this,reactor,speciesStart+reactorStart); 
       //Solving energy equation for preconditioner 
       double dEdt = reactor->evaluateEnergyEquation(t,y,ydot,params);   
        //Solving other variables
        for (unsigned long i = 0; i < speciesStart; i++)
        {
            std::string component = reactor->componentName(i);
            if(!component.compare("temperature"))
            {   
                double dTdt = dEdt/(reactor->mass()*(reactor->getThermoMgr())->cv_mass()); //adjusting energy to get dTdt - K/s
                Cantera::AMP::TemperatureDerivatives(this,reactor,ydot,dTdt,reactorStart+i,speciesStart+reactorStart); //Temperature is index location 2
            }
            else if (!component.compare("mass"))
            {
                Cantera::AMP::NoPrecondition(this,i+reactorStart,i+reactorStart); //setting mass variable element of preconditioner equal to 1
            }
            else if (!component.compare("volume"))
            {
                Cantera::AMP::NoPrecondition(this,i+reactorStart,i+reactorStart); //setting mass variable element of preconditioner equal to 1
            }
            else
            {   std::string errmessage = component+": Adaptive preconditioning is not implemented this variable";
                throw CanteraError("AdaptivePreconditioner::setup", errmessage);
            }
        }
    }

    void AdaptivePreconditioner::initialize(unsigned long nrows,unsigned long ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
        this->matrix.resize(nrows,ncols);
        this->matrix.reserve(nrows*ncols);
    }

    void AdaptivePreconditioner::reset()
    {
        //Do any reset stuff here
        this->matrix.setZero(); //Set all elements to zero
        this->matrix.makeCompressed(); //Compress matrix
    }

    //! Use this function to print and check reactor components
    inline void printReactorComponents(Reactor* reactor)
    {
        for (unsigned long i = 0; i < reactor->neq(); i++)
        {
        std::cout<<reactor->componentName(i)<<std::endl;
        }
    }

    //! This function determines the rate of progress derivatives given a composition of reactants or products
    inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,unsigned long> indexMap, double* omega, double* concentrations, double k_direction, double volume)
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

    //!This function does not precondition the associated equation by assigning it's preconditioner value to a value of 1
    //!@param row the row index of the variable
    //!@param col the column index of the variable
    void NoPrecondition(PreconditionerBase* preconditioner,unsigned long row, unsigned long col)
    {
        preconditioner->setElement(row,col,1); //setting mass variable element of preconditioner equal to 1
    }

    //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
    //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
    //! @param *reactor A pointer to the current reactor being precondition
    //! @param *ydot A pointer to the current data of ydot passed from CVODES
    //! @param meanSpecificHeat The mean specific heat used based on reactor type
    //! @param index The index location of temperature in the state vector
    void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double* ydot, double dTdt, unsigned long index, unsigned long speciesStart)
    {   
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        //Important sizes to the determination of values
        unsigned long numberOfSpecies = kinetics->nTotalSpecies();
        //Perturbation for finite difference of temperature
        double perturbationInverse = 1/std::sqrt(__DBL_EPSILON__);
        //Array pointers for data that is reused
        //net production rates (omega dot)
        double* netProductionRates = new double[numberOfSpecies];
        kinetics->getNetProductionRates(netProductionRates);
        //Adding to preconditioner the species derivatives
        double specTempDerivative;
        for (unsigned long j = 0; j < numberOfSpecies; j++) //column
        {   
            specTempDerivative = (netProductionRates[j]-ydot[j+speciesStart])*perturbationInverse; //+1 because specie index will start after temperature
            preconditioner->setElementByThreshold(index,j+speciesStart,specTempDerivative); //Add by threshold
        }
        //Adding to preconditioner the temperature derivative
        dTdt -= ydot[index];
        dTdt *= perturbationInverse; //Turning dTdt into d(dTdt)/dT
        preconditioner->setElementByThreshold(index,index,dTdt); //Add by threshold
        //Deleting appropriate pointers
        delete[] netProductionRates;
        }

        //! This function determines derivatives of Species with respect to species for jacobian preconditioning;
        //! specifically it determines the derivatives of the rate laws of all species with respect to other species in terms of moles.
        //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
        //! @param *reactor A pointer to the current reactor being precondition
        void SpeciesSpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, unsigned long speciesStart)
        {
        
        preconditioner->setElementByThreshold(0,0,1.0);
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
        double* kForward = new double[numberOfReactions];
        double* kBackward = new double[numberOfReactions];
        //Concentrations of species
        double* concentrations = new double[numberOfSpecies];
        //rate law molar derivatives
        double* rateLawDerivatives = new double[numberOfSpecies*numberOfSpecies];
        //Getting species names
        const std::vector<std::string> *names = &(thermo->speciesNames());
        //Getting species concentrations
        thermo->getConcentrations(concentrations);
        //Getting forward rate constants for calcs
        kinetics->getFwdRateConstants(kForward); 
        //Getting reverse rate constants for calcs
        kinetics->getRevRateConstants(kBackward); 
        //shared_ptr for current reaction in finding Jacobian
        std::shared_ptr<Reaction> currentReaction;
        //Creating map of species indices
        std::map<std::string,unsigned long> indexMap;
        for (unsigned long i = 0; i < numberOfSpecies; i++)
        {
            indexMap[names->at(i)]=i;
        }
        for (unsigned long r = 0; r < numberOfReactions; r++)
        {
            currentReaction=kinetics->getReactionPtr(r);
            //Loop through reactants in current reaction
            reactants = currentReaction->reactants;
            products = currentReaction->products;
            Cantera::AMP::speciesDerivative(reactants,indexMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume());
            Cantera::AMP::speciesDerivative(products,indexMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume()); //Multiply by negative one to change direction
        }

        //Adding to preconditioner
        unsigned long idx;
        for (unsigned long j = 0; j < numberOfSpecies; j++) // column
            {
            for (unsigned long i = 0; i < numberOfSpecies; i++) //row
            {  
            idx = j+i*numberOfSpecies; //Getting flattened index
            preconditioner->setElementByThreshold(i+speciesStart,j+speciesStart,rateLawDerivatives[idx]);//Add by threshold
            }
        }

        //Deleting appropriate pointers
        delete[] kForward;
        delete[] kBackward;
        delete[] concentrations;
        delete[] rateLawDerivatives;
    }
}
