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

#ifndef SENSITIVITYDATASTRUCTURE_HPP_
#define SENSITIVITYDATASTRUCTURE_HPP_

#include <fstream>
#include <vector>
#include <map>
#include "Debug.hpp"

//#include "UblasVectorInclude.hpp"
#include "FileFinder.hpp"


class SensitivityDataStructure
{
protected:

    /**
     * @param fileName  of the data file
     */
    void LoadDataFromFile(std::string fileName);

    std::vector<unsigned> mMask;
    std::vector<unsigned> mExperimentNumbers;
    std::vector<std::vector<double> > mScaleFactors;
    std::map<std::string, unsigned> mMapFromMetadataToColumn;

    void SetUpMap();

public:

    /**
     * Default Constructor
     */
    SensitivityDataStructure(std::string fileName);

    /**
     * Destructor
     */
    ~SensitivityDataStructure();

    /**
     * @return The number of experiments in the data file.
     */
    unsigned GetNumExperiments();

    /**
     * @return The number of parameters in the mask
     */
    unsigned GetMaskLength();

    /**
     * @return  Vector of masks for model parameters
     */
    std::vector<unsigned> GetMask();

    /**
    * @param experimentNumber number of the experiment in this file
    * @return  Number of experiment for entire experimental design
    */
    unsigned GetExperimentNumber(unsigned experimentNumber);

    /**
     * @param experimentNumber number of the experiment in this file
     * @param rParameterName  the cellml metadata name of the parameter to scale (defined in python.pycml/cellml_metadata.py)
     * @return  Vector of scale factors for this experiment
     */
    double GetScaleFactor(unsigned experimentNumber, const std::string& rParameterName);

    /**
     *
     * @param rParameterName
     * @return Whether the parameter has a corresponding column in the data file.
     */
    bool IsParameterInFile(const std::string& rParameterName);

};

#endif // SENSITIVITYDATASTRUCTURE_HPP_
