/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifdef CHASTE_CVODE

#ifndef _TESTSENSITIVITYANALYSISOHARAENDO_HPP_
#define _TESTSENSITIVITYANALYSISOHARAENDO_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>
#include <iostream>

#include <boost/assign/list_of.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include "ohara_rudy_2011Cvode.hpp"
#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"
#include "SensitivityDataStructure.hpp"
#include "CellProperties.hpp"
#include "VectorHelperFunctions.hpp"

class TestSensitivityAnalysisOHaraEndo : public CxxTest::TestSuite
{
public:
    /**
     * This test WILL wipe $CHASTE_TEST_OUTPUT/SensitivityAnalysis/<Model>/
     *
     * Parameters can be defined at the top of this Test
     */
    void TestSensitivityOHaraEndo(void) throw (Exception)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments.
        std::cout << "# " << argc-1 << " arguments supplied.\n" << std::flush;

        if ( !CommandLineArguments::Instance()->OptionExists("--file"))
        {
            std::cerr << "TestSensitivityAnalysis::Please input an argument\n"
                         "* --file  the path of an experimental design file\n"
            			 "* --tag  tag for identifying output (optional)\n"
                         "* --store-data  record voltage and calcium transients (optional)\n";
            return;
        }

        FileFinder exp_des_path(CommandLineArguments::Instance()->GetStringCorrespondingToOption("--file"),
                                 RelativeTo::AbsoluteOrCwd);
        if (!exp_des_path.Exists())
        {
            EXCEPTION("No experiment description file found at :" << exp_des_path.GetAbsolutePath());
        }

        // Get tag
        std::string tag = "";
        if (CommandLineArguments::Instance()->OptionExists("--tag"))
		{
        	tag=CommandLineArguments::Instance()->GetStringCorrespondingToOption("--tag");
		}

        // Set up foldernames for each model and protocol set.
        std::string foldername = "SensitivityAnalysis/OHara2011_endo" +tag;

        // Make above directories. This will NOT wipe an existing folder!
        OutputFileHandler steady_handler(foldername, false);

        // Check to see if action potentials and calcium transients are to be stored
        bool store_data = false;
        if (CommandLineArguments::Instance()->OptionExists("--store-data"))
        {
        	store_data=true;
        }

        // Stimulus is set later, use simple stimulus for now. Solver is not required for a Cvode cell.
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-80, 0.5, 0));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        Cellohara_rudy_2011FromCellMLCvode* p_model = new Cellohara_rudy_2011FromCellMLCvode(p_solver, p_stimulus);


        // The following names are fixed and correspond to metadata tags.
        // We record the default parameter values that the model uses and const them to avoid problems!
        const double default_g_na = p_model->GetParameter("membrane_fast_sodium_current_conductance");
        const double default_g_cal= p_model->GetParameter("membrane_L_type_calcium_current_conductance");
        const double default_g_kr = p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
        const double default_g_ks = p_model->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        const double default_g_k1 = p_model->GetParameter("membrane_inward_rectifier_potassium_current_conductance");
        const double default_j_up = p_model->GetParameter("SR_uptake_current_max");
        const double default_g_to = p_model->GetParameter("membrane_transient_outward_current_conductance");
        const double default_g_naca = p_model->GetParameter("membrane_sodium_calcium_exchanger_current_conductance");

        // Increase the max number of steps Cvode can take.
        p_model->SetMaxSteps(1e8);

        boost::shared_ptr<RegularStimulus> p_default_stimulus =
                             boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction());

        double s1_period = 1500; // This will be overwritten below anyway.
        double s_magnitude = p_default_stimulus->GetMagnitude();
        double s_duration = p_default_stimulus->GetDuration(); // ms
        double s_start = 0;  // ms
        double printing_time_step = 1.0; // ms
        double maximum_time_step; // ms

        // PARAMETERS FOR THE STEADY STATE EXPERIMENTS
        // This is a regular stimulus for steady-state pacing experiments
        boost::shared_ptr<RegularStimulus> p_regular_stimulus(new RegularStimulus(s_magnitude, s_duration, s1_period, s_start));

        if (printing_time_step > s_duration)
        {
            maximum_time_step = s_duration; // excludes case where
        }
        else
        {
            maximum_time_step = printing_time_step;
        }

        boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = p_model->GetSystemInformation();
        unsigned voltage_index = p_ode_info->GetStateVariableIndex("membrane_voltage");
        unsigned calcium_index = p_ode_info->GetStateVariableIndex("cytosolic_calcium_concentration");


        // Load parameter scale factor data
        SensitivityDataStructure scale_factor_data(exp_des_path.GetAbsolutePath());

        // Pacing frequencies of interest.
        std::vector<double> pacing_cycle_lengths = boost::assign::list_of(1500)(1000)(900)(800)(700)(600)(500)(450)(400)(350)(300);

        /**
         * START LOOP OVER EXPERIMENTS
         *
         */
	    std::cout << scale_factor_data.GetNumExperiments() << " experiments in data file.\n";
        std::cout << "Running simulations..." << std::endl << std::flush;

        for(unsigned experiment_idx = 0; experiment_idx < scale_factor_data.GetNumExperiments(); experiment_idx++)
        {
        	// Get the number for this experiment in the overall experimental design
        	unsigned current_experiment = scale_factor_data.GetExperimentNumber(experiment_idx);

        	std::string current_experiment_string = boost::lexical_cast<std::string>(current_experiment);

        	std::string OutputFileName;
        	std::stringstream OutputFileNameStream;
        	OutputFileNameStream << "DYN_" << current_experiment <<".out";
        	OutputFileName = OutputFileNameStream.str();

            // Open output files
            out_stream biomarker_results_file =  steady_handler.OpenOutputFile(OutputFileName);

            // Reset to initial conditions (avoid starting from an unusual state)
            p_model->SetStateVariables(p_model->GetInitialConditions());

            // Set up the parameters for this experiment.
            double gNa_factor = 1;
            double gCaL_factor  = 1;
            double gKr_factor = 1;
            double gKs_factor = 1;
            double gK1_factor = 1;
            double jUp_factor = 1;
            double gTo_factor = 1;
            double gNaCa_factor = 1;

            if (scale_factor_data.IsParameterInFile("membrane_fast_sodium_current_conductance"))
            {
                gNa_factor = scale_factor_data.GetScaleFactor(experiment_idx, "membrane_fast_sodium_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("membrane_L_type_calcium_current_conductance"))
            {
                gCaL_factor = scale_factor_data.GetScaleFactor(experiment_idx, "membrane_L_type_calcium_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("membrane_rapid_delayed_rectifier_potassium_current_conductance"))
            {
                gKr_factor = scale_factor_data.GetScaleFactor(experiment_idx, "membrane_rapid_delayed_rectifier_potassium_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("membrane_slow_delayed_rectifier_potassium_current_conductance"))
            {
                gKs_factor = scale_factor_data.GetScaleFactor(experiment_idx, "membrane_slow_delayed_rectifier_potassium_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("membrane_inward_rectifier_potassium_current_conductance"))
            {
                gK1_factor = scale_factor_data.GetScaleFactor(experiment_idx, "membrane_inward_rectifier_potassium_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("SR_uptake_current_max"))
            {
                jUp_factor = scale_factor_data.GetScaleFactor(experiment_idx,"SR_uptake_current_max");
            }
            if (scale_factor_data.IsParameterInFile("membrane_transient_outward_current_conductance"))
            {
                gTo_factor = scale_factor_data.GetScaleFactor(experiment_idx,"membrane_transient_outward_current_conductance");
            }
            if (scale_factor_data.IsParameterInFile("membrane_sodium_calcium_exchanger_current_conductance"))
            {
                gNaCa_factor = scale_factor_data.GetScaleFactor(experiment_idx,"membrane_sodium_calcium_exchanger_current_conductance");
            }


            std::cout << "gNa factor = " << gNa_factor << "\n" << std::flush;
            std::cout << "gCaL factor = " << gCaL_factor << "\n" << std::flush;
            std::cout << "gKr factor = " << gKr_factor << "\n" << std::flush;
            std::cout << "gKs factor = " << gKs_factor << "\n" << std::flush;
            std::cout << "gK1 factor = " << gK1_factor << "\n" << std::flush;
            std::cout << "jUp factor = " << jUp_factor << "\n" << std::flush;
            std::cout << "gTo factor = " << gTo_factor << "\n" << std::flush;
            std::cout << "gNaCa factor = " << gNaCa_factor << "\n" << std::flush;

            // The following names are fixed and correspond to metadata tags.
            p_model->SetParameter("membrane_fast_sodium_current_conductance",default_g_na*gNa_factor);
            p_model->SetParameter("membrane_L_type_calcium_current_conductance",default_g_cal*gCaL_factor);
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance",default_g_kr*gKr_factor);
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance",default_g_ks*gKs_factor);
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance",default_g_k1*gK1_factor);
            p_model->SetParameter("SR_uptake_current_max",default_j_up*jUp_factor);
            p_model->SetParameter("membrane_transient_outward_current_conductance",default_g_to*gTo_factor);
            p_model->SetParameter("membrane_sodium_calcium_exchanger_current_conductance",default_g_naca*gNaCa_factor);

            // Use optimisations from #1912 - Reset on the next Solve call.
            p_model->SetAutoReset(true);

            /**
             * Start loop over pacing frequencies of interest.
             */
            for (unsigned pacing_idx=0; pacing_idx<pacing_cycle_lengths.size(); ++pacing_idx)
            {
                out_stream data_traces_file;
                if(store_data)
                {
                    std::stringstream file_name_stream;
                    file_name_stream << "DataTraces_" << current_experiment << "_" << pacing_cycle_lengths[pacing_idx]<< "ms.dat";
                    std::string data_traces_file_name = file_name_stream.str();
                    data_traces_file = steady_handler.OpenOutputFile(data_traces_file_name);
                }


            	/**
            	 * STEADY STATE PACING EXPERIMENT
            	 */
                s1_period = pacing_cycle_lengths[pacing_idx];
                p_regular_stimulus->SetPeriod(s1_period);
                p_model->SetStimulusFunction(p_regular_stimulus); // Assign the regular stimulus to the cell's stimulus

                //SteadyStateRunner steady_runner(p_model);
                //steady_runner.RunToSteadyState();
                p_model->Solve(0, 98*s1_period, maximum_time_step);
                // Use optimisations from #1912 - Don't reset on the next solve call.
                p_model->SetAutoReset(false);

                std::vector<double> apd80_vector;
                std::vector<double> apd30_vector;
                std::vector<double> caitd80_vector;
                std::vector<double> caitd30_vector;
                std::vector<double> calcium_peak_vector;
                std::vector<double> delay_vector;

                // Get the solution for a final 2 paces and get more detail....
                unsigned num_repeats = 2u;
                // Prepare somewhere to store the traces
                std::vector<std::vector<double> > store_results(2*num_repeats+1);

                for(unsigned repeats=0; repeats<num_repeats; repeats++)
                {
                    OdeSolution solution = p_model->Solve(0, s1_period, maximum_time_step/10, maximum_time_step/10);

                    // Get voltage and calcium properties
                    std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index); // Voltage should always be zero
                    std::vector<double> cytosolic_calcium = solution.GetVariableAtIndex(calcium_index);

                    if (store_data)
                    {
                        store_results[0] = solution.rGetTimes(); // same on all repeats
                        store_results[2*repeats+1] = voltages;
                        store_results[2*repeats+2] = cytosolic_calcium;
                    }

                    CellProperties voltage_properties(voltages, solution.rGetTimes(),-30);// Threshold for AP is above -70ish.

                    double apd80, apd30, ap_upstroke_time;//, ap_peak;
                    try
                    {
                        apd80 = voltage_properties.GetLastActionPotentialDuration(80);
                        apd30 = voltage_properties.GetLastActionPotentialDuration(30);
                        ap_upstroke_time = voltage_properties.GetTimeAtLastMaxUpstrokeVelocity();
                        //ap_peak = voltage_properties.GetLastPeakPotential();
                    }
                    catch (Exception& e)
                    {
                        if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                e.GetShortMessage()=="No full action potential was recorded")
                        {
                            std::cout << "At Steady pacing of " << s1_period << "ms we did not record any APs\n" << std::flush;
                            apd80 = -1;
                            apd30 = -1;
                            ap_upstroke_time = -1;
                            //ap_peak = -1; // Default error/missing data code for matlab analysis.
                        }
                        else
                        {
                            throw e;
                        }
                    }


                    // Calculate threshold based on max - min of calcium transient
                    double max_ca = -DBL_MAX;
                    double min_ca = DBL_MAX;
                    for (unsigned i=0 ; i<cytosolic_calcium.size(); i++)
                    {
                        if (cytosolic_calcium[i] > max_ca)
                        {
                            max_ca = cytosolic_calcium[i];
                        }
                        if (cytosolic_calcium[i] < min_ca)
                        {
                            min_ca = cytosolic_calcium[i];
                        }
                    }
                    double ca_threshold = min_ca + (max_ca - min_ca)*0.2;

                    CellProperties calcium_properties(cytosolic_calcium, solution.rGetTimes(),ca_threshold);// Threshold for CaiT is above 1e-5 mM

                    double caitd80, caitd30, calcium_upstroke_time, calcium_peak;
                    try
                    {
                        caitd80 = calcium_properties.GetLastActionPotentialDuration(80);
                        caitd30 = calcium_properties.GetLastActionPotentialDuration(30);
                        calcium_upstroke_time = calcium_properties.GetTimeAtLastMaxUpstrokeVelocity();
                        calcium_peak = calcium_properties.GetLastPeakPotential();
                    }
                    catch (Exception& e)
                    {
                        if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                e.GetShortMessage()=="No full action potential was recorded")
                        {
                            std::cout << "At Steady pacing of " << s1_period << "ms we did not record any calcium transients\n" << std::flush;
                            caitd80 = -1;
                            caitd30 = -1;
                            calcium_upstroke_time = -1;
                            calcium_peak = -1; // Default error/missing data code for matlab analysis.
                        }
                        else
                        {
                            throw e;
                        }
                    }
                    double delay;
                    if (calcium_upstroke_time != -1 && ap_upstroke_time!=-1)
                    {
                        delay = calcium_upstroke_time - ap_upstroke_time; // calcium follows upstroke, so this is positive.
                    }
                    else
                    {
                        delay = -1;
                    }

                    apd80_vector.push_back(apd80);
                    apd30_vector.push_back(apd30);
                    caitd80_vector.push_back(caitd80);
                    caitd30_vector.push_back(caitd30);
                    calcium_peak_vector.push_back(calcium_peak);
                    delay_vector.push_back(delay);

                }

                if(store_data)
                {
                    for (unsigned i=0; i<store_results[0].size(); i++)
                    {
                        *data_traces_file << store_results[0][i] << "\t" <<
                                             store_results[1][i] << "\t" <<
                                             store_results[2][i] << "\t" <<
                                             store_results[3][i] << "\t" <<
                                             store_results[4][i] << "\n";
                    }
                    data_traces_file->close();

                }

                // Write out results...
                *biomarker_results_file << pacing_cycle_lengths[pacing_idx] << "\t"  <<
                        apd30_vector[0]  << "\t" << apd30_vector[1] << "\t"  <<
                        apd80_vector[0]  << "\t" << apd80_vector[1] << "\t"  <<
                        caitd30_vector[0]  << "\t" << caitd30_vector[1] << "\t"  <<
                        caitd80_vector[0]  << "\t" << caitd80_vector[1] << "\t"  <<
                        delay_vector[0]  << "\t" << delay_vector[1] <<"\t"  <<
                        calcium_peak_vector[0]  << "\t" << calcium_peak_vector[1] <<"\n" << std::flush;


            }// Pacing Freq
            // Tidy up
            biomarker_results_file->close();

        } // Experiment

    } // Test

};


#endif //_TESTSENSITIVITYANALYSIS_HPP_

#endif //_CHASTE_CVODE

