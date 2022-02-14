/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

// set the tumor cell properties, then call the function 
// to set up the tumor cells 

void create_TH_cells();
void create_CTL_cells();
void create_stroma_cells();
void create_cancer_cells();
void create_cell_types();

void setup_tissue_circle_immune();

void setup_microenvironment( void );  // done

void cancer_cell_proliferation_infection_movement( Cell* pCell, Phenotype& phenotype, double dt ); 

void cell_movement( Cell* pCell, Phenotype& phenotype, double dt );

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt );

void stroma_function( Cell* pCell, Phenotype& phenotype, double dt );

void TH_functions( Cell* pCell, Phenotype& phenotype, double dt);

void CTL_functions( Cell* pCell, Phenotype& phenotype, double dt );

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt );

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant );

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt );

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt );

void attach_cells( Cell* pCell_1, Cell* pCell_2 );

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 );

void CTL_movement( Cell* pCell, Phenotype& phenotype, double dt );

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt );

void virus_induced_lysis( Cell* pCell, Phenotype& phenotype, double dt );

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget );

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell );

std::vector<double> speed_distribution(void);
std::vector<double> go_times_cumulative(void);
std::vector<double> persistance_distribution(void);
std::vector<double> speed_cumulative(void);
	
void immune_cell_placement(void);

