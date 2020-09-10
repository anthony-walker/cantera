#include "cantera/numerics/Preconditioners.h"

namespace Cantera::AMP //Making ASP apart of Cantera namespace
{
    //! Use this function to print and check reactor components
    inline void printReactorComponents(Reactor* reactor)
    {
        for (size_t i = 0; i < reactor->neq(); i++)
        {
        std::cout<<reactor->componentName(i)<<std::endl;
        }
    }

    //! This function determines the rate of progress derivatives given a composition of reactants or products
    inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,size_t> indexMap, double* omega, double* concentrations, double k_direction, double volume)
    { 
        //flattened index for derivatives
        size_t oidx; //index for omega
        size_t sidx; //index for species
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

    //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
    //! @param *preconditioner A pointer to a SparseMatrix Object for preconditioning the system and storing preconditioner values
    //! @param *reactor A pointer to the current reactor being precondition
    //! @param *ydot A pointer to the current data of ydot passed from CVODES
    //! @param meanSpecificHeat The mean specific heat used based on reactor type
    //! @param index The index location of temperature in the state vector
    void TemperatureDerivatives(SparseMatrix *preconditioner,Reactor* reactor, double* ydot, double dTdt, size_t index, size_t speciesStart)
    {
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        //Important sizes to the determination of values
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        //Perturbation for finite difference of temperature
        double perturbationInverse = 1/std::sqrt(UNIT_ROUNDOFF);
        //Array pointers for data that is reused
        //net production rates (omega dot)
        double* netProductionRates = new double[numberOfSpecies];
        kinetics->getNetProductionRates(netProductionRates);
        //Adding to preconditioner the species derivatives
        double specTempDerivative;
        for (size_t j = 0; j < numberOfSpecies; j++) //column
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
        //! @param *preconditioner A pointer to a SparseMatrix Object for preconditioning the system and storing preconditioner values
        //! @param *reactor A pointer to the current reactor being precondition
        void SpeciesSpeciesDerivatives(SparseMatrix *preconditioner,Reactor* reactor, size_t speciesStart)
        {
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        //Getting thermophase object for access to concentrations and species data
        ThermoPhase* thermo=reactor->getThermoMgr();
        //Compositions for reactants and products
        Composition reactants;
        Composition products;
        //Important sizes to the determination of values
        size_t numberOfReactions = kinetics->nReactions();
        size_t numberOfSpecies = kinetics->nTotalSpecies();
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
        std::map<std::string,size_t> indexMap;
        for (size_t i = 0; i < numberOfSpecies; i++)
        {
            indexMap[names->at(i)]=i;
        }
        for (size_t r = 0; r < numberOfReactions; r++)
        {
            currentReaction=kinetics->getReactionPtr(r);
            //Loop through reactants in current reaction
            reactants = currentReaction->reactants;
            products = currentReaction->products;
            Cantera::AMP::speciesDerivative(reactants,indexMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume());
            Cantera::AMP::speciesDerivative(products,indexMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume()); //Multiply by negative one to change direction
        }

        //Adding to preconditioner
        size_t idx;
        for (size_t j = 0; j < numberOfSpecies; j++) // column
            {
            for (size_t i = 0; i < numberOfSpecies; i++) //row
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