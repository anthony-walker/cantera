/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVE_SPARSE_PRECONDITIONER_H
#define ADAPTIVE_SPARSE_PRECONDITIONER_H

const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;
const int PRECONDITIONER_NOT_SET = 0;


//Cantera imports
#include "cantera/numerics/SparseMatrix.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Reaction.h"
//Sundials imports
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif
//Other imports

//Function Library for preconditioniong
namespace Cantera //Making ASP apart of Cantera namespace
{
  namespace AMP //AdaptiveMechanismPreconditioner
  {

    inline void printReactorComponents(Reactor* reactor)
    {
      for (size_t i = 0; i < reactor->neq(); i++)
      {
        std::cout<<reactor->componentName(i)<<std::endl;
      }
    }

    /**
     * 
     * Template Functions
     * 
     **/
    //! This function determines the rate of progress derivatives given a composition of reactants or products
    inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,size_t> indexMap, double* omega, double* concentrations, double k_direction, double volume){ 
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
    /**
     * 
     * Template Functions
     * 
     **/

    //! This function determines derivatives of Species with respect to Temperature for jacobian preconditioning.
    template<class MATTYPE> void SpeciesTemperatureDerivative(SparseMatrix<MATTYPE> *preconditioner,Reactor* reactor, size_t row)
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
        size_t numberOfEquations = reactor->neq();
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        //Array pointers for data that is reused
        double* kForward = new double[numberOfReactions];
        double* kBackward = new double[numberOfReactions];
        //Concentrations of species
        double* concentrations = new double[numberOfEquations];
        //rate law molar derivatives
        double* temperatureDerivatives = new double[numberOfSpecies];
        //Getting species names
        const std::vector<std::string> *names = &(thermo->speciesNames());
        //Getting species concentrations
        thermo->getConcentrations(concentrations);
        //Getting forward rate constants for calcs
        kinetics->getFwdRateConstants(kForward); 
        //Getting reverse rate constants for calcs
        kinetics->getRevRateConstants(kBackward); 
        //shared_ptr for current reaction in finding Jacobian
        ElementaryReaction* currentReaction;
        //rate object to get date for derivative
        Arrhenius rate;
        double dkfdT;
        //Creating map of species indices
        std::map<std::string,size_t> indexMap;
        for (size_t i = 0; i < numberOfSpecies; i++)
        {
          indexMap[names->at(i)]=i;
        }
      
        for (size_t r = 0; r < numberOfReactions; r++)
        {
          currentReaction=(ElementaryReaction*)(kinetics->getReactionPtr(r)).get();
          rate = currentReaction->rate;

          // Loop through reactants in current reaction
          reactants = currentReaction->reactants;
          products = currentReaction->products;
          //Use finite difference here to determine value
          // //Derivative or kf wrt temperature for this reaction
          // dkfdT = kForward[r]/reactor->temperature()*(rate.temperatureExponent()+rate.activationEnergy_R()/(reactor->temperature()*GasConstant));
          // // Forward rate derivative
          // for (std::map<std::string,double>::iterator iterf = reactants.begin(); iterf != reactants.end(); iterf++) //Dependent variable
          // {
            
          // }
        }

        //Adding to preconditioner
        size_t diff = numberOfEquations-numberOfSpecies;
        for (size_t j = 0; j < numberOfSpecies; j++) //column
        {
            preconditioner->setElement(row,j+diff,temperatureDerivatives[j],true);//Add by threshold
        }

        //Deleting appropriate pointers
        delete[] kForward;
        delete[] kBackward;
        delete[] concentrations;
        delete[] temperatureDerivatives;
    }
    //! This function determines derivatives of Species with respect to Pressure for jacobian preconditioning.
    template<class MATTYPE> void SpeciesPressureDerivative(SparseMatrix<MATTYPE> *preconditioner,Reactor* reactor, size_t row)
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
      size_t numberOfEquations = reactor->neq();
      size_t numberOfSpecies = kinetics->nTotalSpecies();
      //Array pointers for data that is reused
      double* kForward = new double[numberOfReactions];
      double* kBackward = new double[numberOfReactions];
      //Concentrations of species
      double* concentrations = new double[numberOfEquations];
      //rate law molar derivatives
      double* volumeDerivatives = new double[numberOfSpecies];
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
      // for (size_t r = 0; r < numberOfReactions; r++)
      // {
      //   currentReaction=kinetics->getReactionPtr(r);
      //   //Loop through reactants in current reaction
      //   reactants = currentReaction->reactants;
      //   products = currentReaction->products;
      //   Cantera::AMP::speciesDerivative(reactants,indexMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume());
      //   Cantera::AMP::speciesDerivative(products,indexMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume()); //Multiply by negative one to change direction
      // }

      //Adding to preconditioner
      size_t diff = numberOfEquations-numberOfSpecies;
      for (size_t j = 0; j < numberOfSpecies; j++) //column
      {
          preconditioner->setElement(row,j+diff,volumeDerivatives[j],true);//Add by threshold
      }

      //Deleting appropriate pointers
      delete[] kForward;
      delete[] kBackward;
      delete[] concentrations;
      delete[] volumeDerivatives;
    }

    //! This function determines derivatives of Species with respect to volume for jacobian preconditioning.
    template<class MATTYPE> void SpeciesVolumeDerivative(SparseMatrix<MATTYPE> *preconditioner,Reactor* reactor, size_t row)
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
      size_t numberOfEquations = reactor->neq();
      size_t numberOfSpecies = kinetics->nTotalSpecies();
      //Array pointers for data that is reused
      double* kForward = new double[numberOfReactions];
      double* kBackward = new double[numberOfReactions];
      //Concentrations of species
      double* concentrations = new double[numberOfEquations];
      //rate law molar derivatives
      double* volumeDerivatives = new double[numberOfSpecies];
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
      // for (size_t r = 0; r < numberOfReactions; r++)
      // {
      //   currentReaction=kinetics->getReactionPtr(r);
      //   //Loop through reactants in current reaction
      //   reactants = currentReaction->reactants;
      //   products = currentReaction->products;
      //   Cantera::AMP::speciesDerivative(reactants,indexMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume());
      //   Cantera::AMP::speciesDerivative(products,indexMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume()); //Multiply by negative one to change direction
      // }

      //Adding to preconditioner
      size_t diff = numberOfEquations-numberOfSpecies;
      for (size_t j = 0; j < numberOfSpecies; j++) //column
      {
          preconditioner->setElement(row,j+diff,volumeDerivatives[j],true);//Add by threshold
      }

      //Deleting appropriate pointers
      delete[] kForward;
      delete[] kBackward;
      delete[] concentrations;
      delete[] volumeDerivatives;
    }

    /**
     * Completed functions
     **/

    //! This function determines derivatives of Species with respect to species for jacobian preconditioning;
    //! specifically it determines the derivatives of the rate laws of all species with respect to other species in terms of moles.
    template<class MATTYPE> void SpeciesSpeciesDerivative(SparseMatrix<MATTYPE> *preconditioner,Reactor* reactor)
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
      size_t numberOfEquations = reactor->neq();
      size_t numberOfSpecies = kinetics->nTotalSpecies();
      //Array pointers for data that is reused
      double* kForward = new double[numberOfReactions];
      double* kBackward = new double[numberOfReactions];
      //Concentrations of species
      double* concentrations = new double[numberOfEquations];
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
      size_t diff = numberOfEquations-numberOfSpecies;
      for (size_t j = 0; j < numberOfSpecies; j++) // column
        {
        for (size_t i = 0; i < numberOfSpecies; i++) //row
        {  
          idx = j+i*numberOfSpecies; //Getting flattened index
          preconditioner->setElement(i+diff,j+diff,rateLawDerivatives[idx],true);//Add by threshold
        }
      }

      //Deleting appropriate pointers
      delete[] kForward;
      delete[] kBackward;
      delete[] concentrations;
      delete[] rateLawDerivatives;
    }
  }
  /*
    End of AdaptiveMechanismPreconditioner namespace
  */
}

#endif