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

#include "SensitivityDataStructure.hpp"
#include "Exception.hpp"

SensitivityDataStructure::SensitivityDataStructure(std::string fileName)
{
	LoadDataFromFile(fileName);

	SetUpMap();
}



SensitivityDataStructure::~SensitivityDataStructure()
{
}

void SensitivityDataStructure::LoadDataFromFile(std::string fileName)
{
   unsigned line_count=1u;
   unsigned numExperimentsInFile;
   unsigned maskLength;
   std::ifstream indata; // indata is like cin
   indata.open(fileName.c_str()); // opens the file
   if(!indata.good())
   { // file couldn't be opened
       EXCEPTION("Couldn't open data file: " + fileName);
   }

   while (indata.good())
   {
       std::string this_line;
       getline(indata, this_line);

       if (this_line=="" || this_line=="\r")
       {
           if (indata.eof())
           {    // If the blank line is the last line carry on OK.
               break;
           }
           else
           {
               EXCEPTION("No data found on this line");
           }
       }
       std::stringstream line(this_line);

       if (line_count==1u)
       {
    	   line >> numExperimentsInFile;
    	   line_count++;
       }
       else if (line_count==2u)
       {
    	   unsigned mask_data;
    	   line >> maskLength;
    	   while(line >> mask_data)
    	   {
    		   mMask.push_back(mask_data);
    	   }
    	   if (maskLength!=mMask.size())
    	   {
    		   EXCEPTION("Length of mask does not equal length of mask vector read");
    	   }
    	   line_count++;

       }
       else
       {
    	   unsigned line_exp_number;
    	   double line_scale_factor;
    	   line >> line_exp_number;
    	   mExperimentNumbers.push_back(line_exp_number);
    	   std::vector<double> line_vector;
    	   while(line >> line_scale_factor)
    	   {
    			   line_vector.push_back(line_scale_factor);
    	   }
    	   if (maskLength!=line_vector.size())
    	   {
    		   EXCEPTION("Length of mask does not equal length of parameter vector read");
    	   }
    	   mScaleFactors.push_back(line_vector);
    	   line_count++;
       }

   }

   if (!indata.eof())
   {
       EXCEPTION("A file reading error occurred");
   }
   if (numExperimentsInFile!=mExperimentNumbers.size() || numExperimentsInFile!=mScaleFactors.size())
   {
	   EXCEPTION("Number of experiments read does not match file header");
   }

}

void SensitivityDataStructure::SetUpMap()
{
    std::map<std::string, unsigned> map_from_metadata_to_mask;
    // Here we define what is hard-coded in Matlab setup - which parameters are what in the mMask.
    map_from_metadata_to_mask["membrane_L_type_calcium_current_conductance"]=0u;
    map_from_metadata_to_mask["membrane_L_type_calcium_current_f_gate_tau"]=1u;
    map_from_metadata_to_mask["membrane_L_type_calcium_current_f2_gate_tau"]=2u;
    map_from_metadata_to_mask["membrane_rapid_delayed_rectifier_potassium_current_conductance"]=3u;
    map_from_metadata_to_mask["membrane_slow_delayed_rectifier_potassium_current_conductance"]=4u;
    map_from_metadata_to_mask["membrane_slow_delayed_rectifier_potassium_current_xs2_gate_tau"]=5u;
    map_from_metadata_to_mask["membrane_inward_rectifier_potassium_current_conductance"]=6u;
    map_from_metadata_to_mask["membrane_sodium_potassium_pump_current_permeability"]=7u;
    map_from_metadata_to_mask["membrane_sodium_calcium_exchanger_current_conductance"]=8u;
    map_from_metadata_to_mask["membrane_transient_outward_current_conductance"]=9u;
    map_from_metadata_to_mask["SR_uptake_current_max"]=10u;

    // Now set up a map from metadata to column index.
    std::map<std::string,unsigned>::iterator it;
    for( it = map_from_metadata_to_mask.begin() ; it != map_from_metadata_to_mask.end(); it++ )
    {
        // If this entry is in the mask add it to the column map
        for( unsigned column=0; column< mMask.size(); column++)
        {
            if (mMask[column]== it->second)
            {
                mMapFromMetadataToColumn[it->first]=column;
            }
        }
    }
}

unsigned SensitivityDataStructure::GetNumExperiments()
{
	return mExperimentNumbers.size();
}

unsigned SensitivityDataStructure::GetMaskLength()
{
	return mMask.size();
}

std::vector<unsigned> SensitivityDataStructure::GetMask()
{
	return mMask;
}

unsigned SensitivityDataStructure::GetExperimentNumber(unsigned experimentNumber)
{
	return mExperimentNumbers[experimentNumber];
}

double SensitivityDataStructure::GetScaleFactor(unsigned experimentNumber, const std::string& rParameterName)
{
    std::map<std::string, unsigned>::iterator it = mMapFromMetadataToColumn.find(rParameterName);
    if (it != mMapFromMetadataToColumn.end())
    {
        return mScaleFactors[experimentNumber][it->second];
    }
    else
    {
        EXCEPTION(rParameterName << " does not correspond to a column of this data file.");
    }
}

bool SensitivityDataStructure::IsParameterInFile(const std::string& rParameterName)
{
    std::map<std::string, unsigned>::iterator it = mMapFromMetadataToColumn.find(rParameterName);
    if (it != mMapFromMetadataToColumn.end())
    {
        return true;
    }
    else
    {
        return false;
    }
}
