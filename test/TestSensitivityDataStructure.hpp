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

#ifndef _TESTSENSITIVITYDATASTRUCTURE_HPP_
#define _TESTSENSITIVITYDATASTRUCTURE_HPP_

#include <cxxtest/TestSuite.h>
//#include <ctime>
//#include <iostream>
#include <vector>

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/shared_ptr.hpp>

#include "Exception.hpp"
#include "Debug.hpp"
//#include "OutputFileHandler.hpp"

#include "SensitivityDataStructure.hpp"
//#include "VectorHelperFunctions.hpp"


class TestSensitivityDataStructure : public CxxTest::TestSuite
{
public:
    /**
     *
     * Parameters can be defined at the top of this Test
     */
    void TestLoadConductanceFile(void) throw (Exception)
    {
    	// Set parameter values for test.
        std::vector<unsigned> mask_values;
        mask_values.push_back(0u);
        mask_values.push_back(3u);
        mask_values.push_back(4u);
        mask_values.push_back(6u);
        mask_values.push_back(8u);
        mask_values.push_back(9u);
        mask_values.push_back(10u);

        std::vector<double> conductances_9576;
        conductances_9576.push_back(0.7090);
        conductances_9576.push_back(1.0229);
        conductances_9576.push_back(0.8163);
        conductances_9576.push_back(0.8684);
        conductances_9576.push_back(0.8395);
        conductances_9576.push_back(0.8527);
        conductances_9576.push_back(0.8393);

        std::vector<double> conductances_9600;
        conductances_9600.push_back(1.0110);
        conductances_9600.push_back(1.2018);
        conductances_9600.push_back(0.7191);
        conductances_9600.push_back(0.7390);
        conductances_9600.push_back(1.1442);
        conductances_9600.push_back(1.2357);
        conductances_9600.push_back(0.8827);

        // Load conductance data
        SensitivityDataStructure sensitivity_data("./projects/JohnW/test/data/Sensitivity/exp_design_uniform_30_384.dat");

        // Check number of experiments
        unsigned num_experiments;
        num_experiments=sensitivity_data.GetNumExperiments();
        TS_ASSERT_EQUALS(num_experiments, 25u);

        // Check length of mask
        unsigned mask_length;
        mask_length = sensitivity_data.GetMaskLength();
        TS_ASSERT_EQUALS(mask_length, 7u);

        // Check mask contents
        std::vector<unsigned> mask;
        mask = sensitivity_data.GetMask();
        TS_ASSERT_EQUALS(mask, mask_values);

        // Check experiment numbers
        unsigned experiment_number1;
        unsigned experiment_number14;
        unsigned experiment_number25;
        experiment_number1 = sensitivity_data.GetExperimentNumber(0u);
        experiment_number14 = sensitivity_data.GetExperimentNumber(13u);
        experiment_number25 = sensitivity_data.GetExperimentNumber(num_experiments-1u);
        TS_ASSERT_EQUALS(experiment_number1, 9576u);
        TS_ASSERT_EQUALS(experiment_number14, 9589u);
        TS_ASSERT_EQUALS(experiment_number25, 9600u);

        // Check data structure contents
        double conductance1_ICaL;
        double conductance25_Jup;

        TS_ASSERT_THROWS_THIS(sensitivity_data.GetScaleFactor(0u, "sausages"),
                              "sausages does not correspond to a column of this data file.");

        TS_ASSERT_EQUALS(sensitivity_data.IsParameterInFile("membrane_L_type_calcium_current_conductance"),true);
        TS_ASSERT_EQUALS(sensitivity_data.IsParameterInFile("sausages"),false);
        TS_ASSERT_EQUALS(sensitivity_data.IsParameterInFile("membrane_L_type_calcium_current_f_gate_tau"),false);

        conductance1_ICaL = sensitivity_data.GetScaleFactor(0u, "membrane_L_type_calcium_current_conductance");
        conductance25_Jup = sensitivity_data.GetScaleFactor(num_experiments-1u, "SR_uptake_current_max");
        TS_ASSERT_DELTA(conductance1_ICaL, conductances_9576[0], 1e-12);
        TS_ASSERT_DELTA(conductance25_Jup, conductances_9600[6], 1e-12);
    }

};


#endif //_TESTSENSITIVITYDATASTRUCTURE_HPP_


