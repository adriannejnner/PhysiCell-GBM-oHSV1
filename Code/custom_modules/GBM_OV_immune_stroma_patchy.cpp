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

#include "./GBM_OV_immune_stroma_patchy.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>


Cell_Definition TH_cell; 
Cell_Definition CTL_cell; 
Cell_Definition cancer_cell; 
Cell_Definition stroma_cell; 

void create_TH_cells( void )
{
	TH_cell = cell_defaults;
	
	TH_cell.name = "TH cell";
	TH_cell.type = 1;
	
	// proliferation 
	TH_cell.functions.cycle_model = Ki67_basic;
	TH_cell.phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	TH_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4"); 	
	TH_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("r_10"); 

	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	TH_cell.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	TH_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	TH_cell.phenotype.motility.migration_speed = parameters.ints("nu");
	
	// cell morphology
	TH_cell.phenotype.geometry.radius = parameters.doubles("R_cell_CD4_CD8");
	TH_cell.phenotype.volume.total = 4/3*3.141*parameters.doubles("R_cell_CD4_CD8")*parameters.doubles("R_cell_CD4_CD8")*parameters.doubles("R_cell_CD4_CD8");;
	TH_cell.phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	TH_cell.phenotype.volume.fluid = TH_cell.phenotype.volume.fluid_fraction*TH_cell.phenotype.volume.total;
	TH_cell.phenotype.volume.solid = TH_cell.phenotype.volume.total-TH_cell.phenotype.volume.fluid;
	TH_cell.phenotype.volume.nuclear = parameters.doubles("V_N_CD4");
	TH_cell.phenotype.volume.nuclear_solid = 23.8;
	TH_cell.phenotype.volume.nuclear_fluid = TH_cell.phenotype.volume.nuclear - TH_cell.phenotype.volume.nuclear_solid;
	TH_cell.phenotype.volume.cytoplasmic = TH_cell.phenotype.volume.total - TH_cell.phenotype.volume.nuclear;
	TH_cell.phenotype.volume.cytoplasmic_fluid = TH_cell.phenotype.volume.fluid_fraction*TH_cell.phenotype.volume.cytoplasmic;
	TH_cell.phenotype.volume.cytoplasmic_solid = TH_cell.phenotype.volume.cytoplasmic-TH_cell.phenotype.volume.cytoplasmic_fluid;
	TH_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.05;
	TH_cell.phenotype.volume.target_solid_cytoplasmic = TH_cell.phenotype.volume.cytoplasmic_solid;
	TH_cell.phenotype.volume.target_solid_nuclear = TH_cell.phenotype.volume.nuclear_solid;
	TH_cell.phenotype.volume.target_fluid_fraction = TH_cell.phenotype.volume.fluid_fraction;
	TH_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = TH_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update function
	TH_cell.functions.update_phenotype = TH_functions;
	
	return;
}
void create_CTL_cells( void )
{
	CTL_cell = cell_defaults;
	
	CTL_cell.name = "CTL cell";
	CTL_cell.type = 3;
	
	// proliferation 
	CTL_cell.functions.cycle_model = Ki67_basic;
	CTL_cell.phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive); 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");//7.2206*1e-5; 
	CTL_cell.phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("r_10");//0.00143; 	
	
	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	CTL_cell.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	CTL_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	CTL_cell.phenotype.motility.migration_speed = parameters.doubles("nu");
	
	//cell morphology
	CTL_cell.phenotype.geometry.radius = parameters.doubles("R_cell_CD4_CD8");
	CTL_cell.phenotype.volume.total = 4/3*3.141*parameters.doubles("R_cell_CD4_CD8")*parameters.doubles("R_cell_CD4_CD8")*parameters.doubles("R_cell_CD4_CD8");
	CTL_cell.phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	CTL_cell.phenotype.volume.fluid = CTL_cell.phenotype.volume.fluid_fraction*CTL_cell.phenotype.volume.total;
	CTL_cell.phenotype.volume.solid = CTL_cell.phenotype.volume.total-CTL_cell.phenotype.volume.fluid;
	CTL_cell.phenotype.volume.nuclear = parameters.doubles("V_N_CD8");
	CTL_cell.phenotype.volume.nuclear_solid = 24.1;
	CTL_cell.phenotype.volume.nuclear_fluid = CTL_cell.phenotype.volume.nuclear - CTL_cell.phenotype.volume.nuclear_solid;
	CTL_cell.phenotype.volume.cytoplasmic = CTL_cell.phenotype.volume.total - CTL_cell.phenotype.volume.nuclear;
	CTL_cell.phenotype.volume.cytoplasmic_fluid = CTL_cell.phenotype.volume.fluid_fraction*CTL_cell.phenotype.volume.cytoplasmic;
	CTL_cell.phenotype.volume.cytoplasmic_solid = CTL_cell.phenotype.volume.cytoplasmic-CTL_cell.phenotype.volume.cytoplasmic_fluid;
	CTL_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.03;
	CTL_cell.phenotype.volume.target_solid_cytoplasmic = CTL_cell.phenotype.volume.cytoplasmic_solid;
	CTL_cell.phenotype.volume.target_solid_nuclear = CTL_cell.phenotype.volume.nuclear_solid;
	CTL_cell.phenotype.volume.target_fluid_fraction = CTL_cell.phenotype.volume.fluid_fraction;
	CTL_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = CTL_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// cell update phenotype		
	CTL_cell.functions.update_phenotype = CTL_functions;
	
	
	return;
}void create_stroma_cells( void )
{
	stroma_cell = cell_defaults;
	
	stroma_cell.name = "stroma cell";
	stroma_cell.type = 4;
	
	// turn off proliferation 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live); 
	stroma_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0; 

	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	stroma_cell.phenotype.secretion.uptake_rates[virus_index] = parameters.doubles("u_s");
	stroma_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	stroma_cell.phenotype.motility.migration_speed = 0;
	stroma_cell.phenotype.motility.is_motile = false;
	stroma_cell.phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;
	
	// cell morphology
	stroma_cell.phenotype.geometry.radius = parameters.doubles("R_cell_stroma");//7.5;
	stroma_cell.phenotype.volume.total = 4/3*3.141*parameters.doubles("R_cell_stroma")*parameters.doubles("R_cell_stroma")*parameters.doubles("R_cell_stroma");//1767;
	stroma_cell.phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	stroma_cell.phenotype.volume.fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.total;
	stroma_cell.phenotype.volume.solid = stroma_cell.phenotype.volume.total-stroma_cell.phenotype.volume.fluid;
	stroma_cell.phenotype.volume.nuclear = parameters.doubles("V_N_stroma");
	stroma_cell.phenotype.volume.nuclear_solid = 125;
	stroma_cell.phenotype.volume.nuclear_fluid = stroma_cell.phenotype.volume.nuclear - stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.cytoplasmic = stroma_cell.phenotype.volume.total - stroma_cell.phenotype.volume.nuclear;
	stroma_cell.phenotype.volume.cytoplasmic_fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.cytoplasmic;
	stroma_cell.phenotype.volume.cytoplasmic_solid = stroma_cell.phenotype.volume.cytoplasmic-stroma_cell.phenotype.volume.cytoplasmic_fluid;
	stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 2.53;
	stroma_cell.phenotype.volume.target_solid_cytoplasmic = stroma_cell.phenotype.volume.cytoplasmic_solid;
	stroma_cell.phenotype.volume.target_solid_nuclear = stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.target_fluid_fraction = stroma_cell.phenotype.volume.fluid_fraction;
	stroma_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update phenotype
	stroma_cell.functions.update_phenotype = stroma_function;
	
	return;
}
void create_cancer_cells( void )
{
	
	cancer_cell = cell_defaults;
	
	//setting cycle model to live
	cancer_cell.phenotype.cycle.sync_to_cycle_model( live ); 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	cancer_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("beta");
	
	// cell actions
	static int virus_index = microenvironment.find_density_index( "virus");
	cancer_cell.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cancer_cell.phenotype.secretion.secretion_rates[virus_index] = 0.0;
	cancer_cell.phenotype.motility.migration_speed = 0.0;
	
	// cell morphology
	cancer_cell.phenotype.geometry.radius = parameters.doubles("R_cell_GBM");
	cancer_cell.phenotype.volume.total = 4/3*3.1416*parameters.doubles("R_cell_GBM")*parameters.doubles("R_cell_GBM")*parameters.doubles("R_cell_GBM");
	cancer_cell.phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	cancer_cell.phenotype.volume.fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.total;
	cancer_cell.phenotype.volume.solid = cancer_cell.phenotype.volume.total-cancer_cell.phenotype.volume.fluid;
	cancer_cell.phenotype.volume.nuclear = parameters.doubles("V_N_GBM");
	cancer_cell.phenotype.volume.nuclear_solid = 185;
	cancer_cell.phenotype.volume.nuclear_fluid = cancer_cell.phenotype.volume.nuclear - cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.cytoplasmic = cancer_cell.phenotype.volume.total - cancer_cell.phenotype.volume.nuclear;
	cancer_cell.phenotype.volume.cytoplasmic_fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.cytoplasmic;
	cancer_cell.phenotype.volume.cytoplasmic_solid = cancer_cell.phenotype.volume.cytoplasmic-cancer_cell.phenotype.volume.cytoplasmic_fluid;
	cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic/cancer_cell.phenotype.volume.nuclear;//6.0321;
	cancer_cell.phenotype.volume.target_solid_cytoplasmic = cancer_cell.phenotype.volume.cytoplasmic_solid;
	cancer_cell.phenotype.volume.target_solid_nuclear = cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.target_fluid_fraction = cancer_cell.phenotype.volume.fluid_fraction;
	cancer_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	cancer_cell.phenotype.volume.calcified_fraction = 0; 
	cancer_cell.phenotype.volume.calcification_rate = 0;
	cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.35*cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength;
	
	//update phenoypt
	cancer_cell.functions.update_phenotype = cancer_cell_proliferation_infection_movement;
	
	cancer_cell.name = "cancer cell";
	cancer_cell.type = 2; 
	
	return;
}

void create_cell_types( void )
{	
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	//cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	//setting cycle model to live with no proliferation as default
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;
	
 	// turn off death model 
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 

	// no cell velocity in default
	cell_defaults.phenotype.motility.migration_speed = 0.0;//0.05;
		
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 ); // amount of intracellular virus
	cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 ); // how long cells will persist in move or stop phenotype
	cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached 
	
	Parameter<double> paramD;
	paramD = parameters.doubles[ "elastic_coefficient" ]; 
	cell_defaults.custom_data.add_variable( "elastic coefficient" , paramD.units, paramD.value ); 
		
	static int virus_index = microenvironment.find_density_index( "virus");
	static int wall_index = microenvironment.find_density_index( "wall");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = parameters.doubles("rhostar_virus");//;//10; 
	cell_defaults.phenotype.molecular.fraction_released_at_death[virus_index] = 1;//1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[wall_index] = 1;//1;
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
	cell_defaults.phenotype.molecular.fraction_released_at_death[chemokine_index] = 1;//1;
	cell_defaults.phenotype.secretion.secretion_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[chemokine_index] = parameters.doubles("rhostar_chemokine");//5; 
	
	// update cell and phenotype based on virus dynamics only
		cell_defaults.functions.update_phenotype = cancer_cell_proliferation_infection_movement;
	cell_defaults.phenotype.motility.is_motile = true; 
		
	cell_defaults.name = "holder cell"; 
	cell_defaults.type = 0; 
	
	//create cell types
	create_cancer_cells();
	create_CTL_cells();
	create_TH_cells();
	create_stroma_cells();
		
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters
	
	// make sure not override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	initialize_microenvironment(); 	

	static int virus_index = microenvironment.find_density_index( "virus");
	static int wall_index = microenvironment.find_density_index( "wall");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")+10)*(parameters.doubles("R")+10))
		{	
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 1;
			microenvironment(n)[chemokine_index] = 0;
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")-10)*(parameters.doubles("R")-10))
		{
			microenvironment(n)[virus_index] = 0;//parameters.doubles("initial_virus_density");
			microenvironment(n)[wall_index] = 3.5;
			microenvironment(n)[chemokine_index] = 0;
		}
		else if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>250*250 )
		{	
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 5;
			microenvironment(n)[chemokine_index] = 0;
				
		}
		else
		{   
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 10;
			microenvironment(n)[chemokine_index] = 0;
		}
	}
	
	// locations of dense patches
	double x_offset_patch_1 = 173.38;
	double y_offset_patch_1 = -270.5;
	double x_offset_patch_2 = 41.32;
	double y_offset_patch_2 = 1015.1;
	double x_offset_patch_3 = -757.84;
	double y_offset_patch_3 = 128.4;
	double x_offset_patch_4 = 899.79;
	double y_offset_patch_4 = -391.6;
	double x_offset_patch_5 = -9.99;
	double y_offset_patch_5 = -859.9;

	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		if( (x_offset_patch_1-ECMdense[0])*(x_offset_patch_1-ECMdense[0])+(y_offset_patch_1-ECMdense[1])*(y_offset_patch_1-ECMdense[1])<20*20)// Centre of patch 1
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_2-ECMdense[0])*(x_offset_patch_2-ECMdense[0])+(y_offset_patch_2-ECMdense[1])*(y_offset_patch_2-ECMdense[1])<20*20)// Centre of patch 2
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_3-ECMdense[0])*(x_offset_patch_3-ECMdense[0])+(y_offset_patch_3-ECMdense[1])*(y_offset_patch_3-ECMdense[1])<20*20) //Centre of patch 3
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_4-ECMdense[0])*(x_offset_patch_4-ECMdense[0])+(y_offset_patch_4-ECMdense[1])*(y_offset_patch_4-ECMdense[1])<20*20) // Centre of patch 4
		{	
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_5-ECMdense[0])*(x_offset_patch_5-ECMdense[0])+(y_offset_patch_5-ECMdense[1])*(y_offset_patch_5-ECMdense[1])<20*20) // Centre of patch 5
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		
	}
	return; 
}

void setup_tissue_circle_immune( void )
{
	double Radius = parameters.doubles("R");
		
	std::vector<double> go_times_cumul(8);
	std::vector<double> persistence_times_vec(12);
	std::vector<double> speed_cumul(12);
	std::vector<double> speed_vec(12);
	go_times_cumul = go_times_cumulative();
	persistence_times_vec = persistance_distribution();
	speed_cumul = speed_cumulative();
	speed_vec = speed_distribution();
		
	immune_cell_placement();
	
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
		
	double GBM_NO_dense = parameters.ints("N_GBM_dense");
	double stroma_NO_dense = parameters.ints("N_stroma_dense");
	double GBM_NO_sparse = parameters.ints("N_GBM_sparse");
	double stroma_NO_sparse = parameters.ints("N_stroma_sparse"); 
			
	double fragment_radius = 1267.7;
	double area_of_fragment = 3.141*fragment_radius*fragment_radius;
	
	double patch_radius_1 = 317.93*parameters.doubles("kappa");
	double x_offset_patch_1 = 173.38;
	double y_offset_patch_1 = -270.5;
	double area_of_patch_1 = 3.141*patch_radius_1*patch_radius_1;
	double no_GBM_cells_patch_1 = GBM_NO_dense/area_of_fragment*area_of_patch_1;
	
	double patch_radius_2 = 638.2*parameters.doubles("kappa");
	double x_offset_patch_2 = 41.32;
	double y_offset_patch_2 = 1015.1;
	double area_of_patch_2 = 3.141*patch_radius_2*patch_radius_2;
	double no_GBM_cells_patch_2 = GBM_NO_dense/area_of_fragment*area_of_patch_2;
	
	double patch_radius_3 = (407.26+42.5)*parameters.doubles("kappa");
	double x_offset_patch_3 = -757.84;
	double y_offset_patch_3 = 128.4;
	double area_of_patch_3 = 3.141*patch_radius_3*patch_radius_3;
	double no_GBM_cells_patch_3 = GBM_NO_dense/area_of_fragment*area_of_patch_3;
	
	double patch_radius_4 = 307*parameters.doubles("kappa");
	double x_offset_patch_4 = 899.79;
	double y_offset_patch_4 = -391.6;
	double area_of_patch_4 = 3.141*patch_radius_4*patch_radius_4;
	double no_GBM_cells_patch_4 = GBM_NO_dense/area_of_fragment*area_of_patch_4;
	
	double patch_radius_5 = 407.67*parameters.doubles("kappa");
	double x_offset_patch_5 = -9.99;
	double y_offset_patch_5 = -859.9;
	double area_of_patch_5 = 3.141*patch_radius_5*patch_radius_5;
	double no_GBM_cells_patch_5 = GBM_NO_dense/area_of_fragment*area_of_patch_5;
	
	double no_GBM_cells_remaining = GBM_NO_sparse*(1-parameters.doubles("proportion"));

	//GBM cells 	
	for( int i=0; i<no_GBM_cells_remaining; i++ )
	{	
		double r_cell = sqrt(UniformRandom())*patch_radius_1;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp)+x_offset_patch_1;
		y = r_cell*sin(alp)+y_offset_patch_1;
		
		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))<patch_radius_1 ||
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))<patch_radius_2 ||
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))<patch_radius_3 ||
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))<patch_radius_4 ||
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))<patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		
		pCell = create_cell( cancer_cell );		
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_g");//distribution2(generator);//distribution2(generator);
		
		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
		
		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();
			
			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pCell->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		} 
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP
	
		double go_stop_var = UniformRandom();
		
		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}
			
		}
	}
	
	double no_stroma_cells_remaining = stroma_NO_sparse*(1-parameters.doubles("proportion"));

	for( int l=0; l<no_stroma_cells_remaining; l++ )
	{	
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);
		
		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))<patch_radius_1 ||
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))<patch_radius_2 ||
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))<patch_radius_3 ||
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))<patch_radius_4 ||
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))<patch_radius_5 ) //cell is located within the patch keep trying to get new position   
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		pCell = create_cell( stroma_cell );
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_s");//distribution2(generator);//distribution2(generator);
		
	}
	

	double no_GBM_cells_inpatch = GBM_NO_dense*(parameters.doubles("proportion"));

	//GBM cells 	
	for( int i=0; i<no_GBM_cells_inpatch; i++ )
	{	
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);
		
		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))>patch_radius_1 &&
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))>patch_radius_2 &&
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))>patch_radius_3 &&
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))>patch_radius_4 &&
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))>patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		
		pCell = create_cell( cancer_cell );		
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_g");//distribution2(generator);//distribution2(generator);
				
		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
		
		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();
			
			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pCell->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		} 
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP
	
		double go_stop_var = UniformRandom();
		
		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}
			
		}
	}
	
	double no_stroma_cells_inpatch = stroma_NO_dense*(parameters.doubles("proportion"));

	for( int l=0; l<no_stroma_cells_inpatch; l++ )
	{	
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);
		
		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))>patch_radius_1 &&
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))>patch_radius_2 &&
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))>patch_radius_3 &&
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))>patch_radius_4 &&
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))>patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		pCell = create_cell( stroma_cell );
		pCell->assign_position( x , y , 0.0 );
		
		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_s");//distribution2(generator);//distribution2(generator);
				
	}
	
	return;
	
}

void cancer_cell_proliferation_infection_movement( Cell* pCell, Phenotype& phenotype, double dt )
{

	double R = parameters.doubles("R_cell_GBM");
	double SA = 4*3.1416*R*R;
	double s = parameters.doubles("xi");
	double pressure = 6*(1-1/(2*R)*s)*(1-1/(2*R)*s);
	double pressure_scale = 0.027288820670331; 
	double max_pressure = pressure/pressure_scale;
	
	//tumour cell proliferation
	if(pCell->type ==2)
	{	
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("beta");

		if( pCell->state.simple_pressure*pCell->state.simple_pressure>max_pressure) // if cell under too much pressure -> no proliferation
		{pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;}
	
	}
	
	cell_movement( pCell, phenotype, dt);
	infection_dynamics( pCell, phenotype, dt );
		
	return;
	
}

void cell_movement( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	//cell movement
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
	
	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; // 1 = go, 2 = stop
	
	std::vector<double> go_times_cumul(8);
	std::vector<double> persistence_times_vec(12);
	std::vector<double> speed_cumul(12);
	std::vector<double> speed_vec(12);
	go_times_cumul = go_times_cumulative();
	persistence_times_vec = persistance_distribution();
	speed_cumul = speed_cumulative();
	speed_vec = speed_distribution();
	
	if( wall_amount<2 & pCell->type == 2 )
	{
		pCell->phenotype.motility.migration_speed = 0;
	}
	else if(pCell->type != 2)
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
	}
	else		
	{
		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
		{
			// assign new type (stop = 2, or go = 1)
			double new_type_rand = UniformRandom();
			if(new_type_rand<=0.5)// GO
			{
				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type
				
				double speed_var = UniformRandom();
			
				for( int k=0; k<12; )
				{
					if( speed_var> speed_cumul[k] )
					{k++;}
					else
					{
						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
						k = 12;
					}
				}
			} 
			else
			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type

			// assign persistence time - needs to be a real time!
			double go_stop_var = UniformRandom();
			for( int j=0; j<8; )
			{
				if( go_stop_var> go_times_cumul[j] )
				{j++;}
					else
				{
					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
					j = 8;
				}
			}
		}
			
	}
	return;
}

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
			
	//cell infection
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
		
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double u = parameters.doubles("u_g");
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double nstar = parameters.doubles("m_half");//10;//infection threshol
	double alp = parameters.doubles("alpha");//1000;//virus burst number
	double nu = parameters.doubles("gamma");
	double pmax = parameters.doubles("rho_max");//0.0125;
	
	if( n>1)
	{
	//	std::cout<<"m_i = "<<n<<" rho_virus = "<<p<<" u_g = "<<u<<" mstar = "<<nstar<<" alph = "<<alp<<" gamma = "<<nu<<" rho_max = "<<pmax<<std::endl;
	}
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];		
		
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
		if(p<pmax)
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(n/Vvoxel+nstar/Vvoxel);}
		else
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*pmax*pmax/(n/Vvoxel+nstar/Vvoxel)/p;}
		
		if( n > 1 && n <= alp) // update amount inside due to replication 
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n+dt*(nu*n); 
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		
		}
		else if( n > alp-1)
		{
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
			pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->start_death( apoptosis_model_index );	
		}
		else if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<0 )
		{std::cout<<"Negative intracellular virus!! m_i = "<<pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<<std::endl;}
	}
	if( pCell->phenotype.death.dead == true && pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index]>0)
	{
		virus_induced_lysis(pCell, phenotype, dt );
	}
	
	return;
}

void virus_induced_lysis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double pstar = pCell->phenotype.secretion.saturation_densities[virus_signal_index];
	double delta_V = parameters.doubles("delta_V");//0.1466;
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	
	if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]> 1 ) 
	{
			double amount_to_add = (n-n*exp(-delta_V*dt))/Vvoxel;
			if( amount_to_add > pstar-p )
			{
				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
			}
			else
			{
				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n*exp(-delta_V*dt);
			}	
	}

	return;
}

void stroma_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	if(pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]>1e5)
	{
		pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
		pCell->custom_data.variables[intracellular_virus_index].value = 0;
	}
	else
	{
		pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	}
	
	return;
}


void TH_functions( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	std::vector<double> ae_ini(3);
    
	if( wall_amount<2 )// Make TH cells that have left the diameter of the tumour turn around
	{
		
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
		ae_ini = -1*pCell->position;
		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
		
	}
	else 
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
		pCell->phenotype.motility.migration_bias = 0;
	}
	
	// TH secretion of cytokines
	std::vector<Cell*> nearby = pCell->cells_in_my_container();
	
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	Cell* pC = NULL;
	bool stop = false;
	int i=0;
	
	double nstar = parameters.doubles("m_half");

	while( !stop && i < nearby.size() )
	{
		pC = nearby[i];
		if( pC->phenotype.molecular.internalized_total_substrates[virus_signal_index] > nstar &&
			pC->phenotype.death.dead == false &&
			pC != pCell && pC->type != 4)
			{ stop = true; }
		
		i++;
	
		if( stop == false )
			{ pC = NULL; }
	}
	
	if( pC )
	{
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = parameters.doubles("S_chemokine_CD4");//0.0417;//0.1;	
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4")*parameters.doubles("TH_prolif_increase_due_to_stimulus"); 	
		pCell->phenotype.motility.migration_speed = 0.1;
	}
	else
	{
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;//0.1;	
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4"); 	
	}
	return;
}

void CTL_functions( Cell* pCell, Phenotype& phenotype, double dt )
{	

	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment lifetime" ); 
	
	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
		bool dettach_me = false; 
		
		if( immune_cell_attempt_apoptosis( pCell, pCell->state.neighbors[0], dt ) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] ); 
			dettach_me = true; 
		}
			
		if( dettach_me )
		{
			
			dettach_cells( pCell, pCell->state.neighbors[0] ); 
			phenotype.motility.is_motile = true; 
			CTL_movement( pCell, phenotype, dt);
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8")*parameters.doubles("CTL_prolif_increase_due_to_stimulus"); 	
		}
		return; 
	}
	
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		phenotype.motility.is_motile = false; 
		return; 
	}
	
	phenotype.motility.is_motile = true; 
	CTL_movement( pCell, phenotype, dt);
	
	return;
}
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	static double max_attachment_distance = parameters.doubles("d_attach"); // 18.0; 
	
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double internal_virus = pTarget->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	
	double kill_time = parameters.doubles("tau"); // how long the cell needs to attach for  the infected cell to be killed
	double nstar = parameters.doubles("m_half"); // how long the cell needs to attach for  the infected cell to be killed
	
	if( internal_virus > nstar && pTarget->phenotype.death.dead == false && pTarget->type!=4)
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 
	
		attach_cells( pAttacker, pTarget ); 
		pAttacker->custom_data[attach_lifetime_i] = PhysiCell_globals.current_time + kill_time;
		
		return true; 
	}
	
	return false; 
}

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return; 
}

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

void CTL_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	
	std::vector<double> ae_ini(3);
    	
	// TH movement
	if( wall_amount<2 )
	{
		std::cout<<"outside domain "<<pCell->phenotype.motility.migration_speed<<std::endl;
		pCell->phenotype.motility.migration_bias = 1;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(wall_index);
		ae_ini = -1*pCell->position;
	
		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
		
		return;
	}
	else if(chemokine_amount>1e-8)// sample chemotaxis gradient and random walk in that direction
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu")+(parameters.doubles("nu_max")-parameters.doubles("nu"))*(chemokine_amount/(parameters.doubles("nu_star")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("b_CD8");//0.85;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");
		
		return;
	}
	else if(chemokine_amount>1e-3)
	{
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8")*parameters.doubles("CTL_prolif_increase_due_to_stimulus");//7.9026*1e-5*1e1; 	
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu")+(parameters.doubles("nu_max")-parameters.doubles("nu"))*(chemokine_amount/(parameters.doubles("nu_star")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("b_CD8");
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
		
		return;
	}
	else
	{
		pCell->phenotype.motility.migration_bias = 0;
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");
		
		return;
	}
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	
	// CTL kills cell if it has been attached for enough time
	if( pAttacker->custom_data[attach_lifetime_i] < PhysiCell_globals.current_time )
	{ 
		return true; 
	}
	return false; 
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	static int virus_index = microenvironment.find_density_index( "virus" ); 
	
	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	pTarget->phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;
	return true; 
}

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int v_index = microenvironment.find_density_index( "virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	double p_min = 1;
	double p_max = parameters.doubles("alpha");
	
	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	if(n_I>10000 && pCell->type ==2)
	{std::cout<<" High viral infection m_i = "<<n_I<<std::endl;}
	
	if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			
		if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{ 
					output[0] = "orange";
					output[1] = "orange";
					output[2] = "coral";
					output[3] = "coral";
					return output; 
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	
					output[0] = "darkred";
					output[1] = "darkred";
					output[2] = "firebrick";
					output[3] = "firebrick";
			}	
		}
		if( pCell->type == 3)
		{
			if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{   
			    output[0] = "aquamarine";
				output[1] = "lightsteelblue";
				output[2] = "lightskyblue";
				output[3] = "aquamarine";
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{	
		
				output[0] = "darkslateblue";
				output[1] = "darkblue";
				output[2] = "darkblue";
				output[3] = "aquamarine";
			}
		}
		if( pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if( n_I>1)
			{
				double p_min = 1;
				int oncoprotein1 = (int) round((210-139)*(1.0/(p_max-p_min)) * (n_I-1)); 
				int oncoprotein2 = (int) round((180-69)*(1.0/(p_max-p_min)) * (n_I-1)); 
				int oncoprotein3 = (int) round((140-19)*(1.0/(p_max-p_min)) * (n_I-1)); 
				
				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", 210-oncoprotein1, 180-oncoprotein2, 140-oncoprotein3); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
				
				output[0].assign( szTempString );
				output[1]="brown";
				output[2].assign( szTempString );
				output[3]="brown";
				
				return output;
			}
			else
			{
					output[0] = "rgb(104, 55, 99)";//"orchid";//"rgb(255,230,230)";
					output[1] = "rgb(104, 55, 99)";
					output[2] = "rgb(85, 50, 70)";//"plum";//"rgb(255,230,230)";
					output[3] = "rgb(85, 50, 70)";

					return output; 
			}
		}
	if( pCell->phenotype.death.dead == true )
	{ 
			output[0] = "rgb(255, 255, 224)";
			output[1] = "rgb(255, 255, 224)";
			output[2] = "rgb(255, 228, 181)";
			output[3] = "rgb(255, 228, 181)";
		
		return output; 
	}
	if( pCell->type == 4)
	{
		output[0] = "rgb(234, 172, 199)";
		output[1] = "rgb(234, 172, 199)";
		output[2] = "rgb(243, 186, 211)";
		output[3] = "rgb(243, 186, 211)";
	}
	return output;
}

std::vector<double> go_times_cumulative()
{
	// setting up distributions for movement and persistance of cells
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	return go_times_cumul;
}

std::vector<double> persistance_distribution()
{	
	std::vector<double> persistance_times_vec(8);
    persistance_times_vec[0] = 0;
	persistance_times_vec[1] = 30;
	persistance_times_vec[2] = 60;
	persistance_times_vec[3] = 90;
	persistance_times_vec[4] = 120;
	persistance_times_vec[5] = 150;
	persistance_times_vec[6] = 180;
	persistance_times_vec[7] = 240;
	
	return  persistance_times_vec;
}

std::vector<double> speed_cumulative()
{	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;

	return speed_cumul;
}
std::vector<double> speed_distribution()	
{	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
	
	return speed_vec;
}

void immune_cell_placement(void)
{
	double x = 0.0;
	double y = 0.0;
	
	Cell* pCell = NULL; 
	
	//UNPROLIFERATIVE TH CELLS
	x = -18.0725;
	y = 922.2948;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 954.1439;
	y = 56.7601;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	// /
    x = -115.9499;
	y = -349.9367;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
    x = 773.5733;
	y = -613.7738;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
    x = 753.4041;
	y = -219.3846;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///
	x = -292.5135;
	y = -91.2428;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
	x = 449.1283;
	y = 407.6992;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 452.4130;
	y = 616.9196;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  60.4108;
	y =  -1021.8;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 534.5535;
	y =  -882.0731;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  15.5367;
	y = 400.0095;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 412.5616;
	y = -921.4520;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  64.9108;
	y = 1021.7;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 669.1877;
	y = -318.7248;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -483.3989;
	y = 759.6808;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 134.3949;
	y = 347.2437;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -47.5182;
	y = 611.5124;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -738.2632;
	y = -659.5292;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -937.1317;
	y = 156.6883;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -614.4497;
	y = 789.3966;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  437.2258;
	y = -793.7856;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -170.1507;
	y = -132.5528;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -885.0803;
	y = -277.3879;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  838.1983;
	y = -501.9633;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -260.7010;
	y =  604.0757;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  103.6960;
	y = -899.8739;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -23.0385;
	y = -922.1671;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -458.5337;
	y = 474.7504;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -765.6844;
	y = -326.6257;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  360.4622;
	y = 228.2000;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  771.4487;
	y = 319.7186;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -195.5453;
	y = -5.9414;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  114.0970;
	y = -561.2747;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///
	x =  203.4536;
	y = -89.5478;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  235.8782;
	y = 265.3679;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -868.0050;
	y = -405.4009;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
		
	x = -811.6261;
	y = 187.1180;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	//PROLIFERATIVE TH CELLS
	x = 454.1934;
	y = -46.7240;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	
	x = -427.8415;
	y = 616.7429;
	pCell = create_cell( TH_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	
	//UNPROLIFERATIVE CTL CELLS
	x = 27.9074;
	y = 119.1812;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 225.7304;
	y = -728.5426;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = -258.8621;
	y =  265.2323;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -809.2029;
	y =  -60.6306;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 473.7135;
	y =  887.8817;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = -366.2289;
	y = -332.2614;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -197.1932;
	y =  446.5596;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -393.1067;
	y =  -597.4458;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -264.8455;
	y =  -840.0492;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x =  -62.2700;
	y = -690.4383;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 968.1414;
	y =  -363.0428;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 592.2715;
	y =   125.8000;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 231.3138;
	y =  809.8071;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 547.2984;
	y =  -545.0700;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 215.9077;
	y =  557.6151;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = 141.7866;
	y = -321.7799;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
	x = -557.2198;
	y = -155.5629;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 809.8675;
	y =  128.2654; 
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 682.2515;
	y =  547.2600;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = -717.2219;
	y =  432.8312;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	x = 588.0764;
	y =  345.4055;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	
	// PROLIFERATIVE CTL cell
	x = 736.2712;
	y = -90.0571;
	pCell = create_cell( CTL_cell );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;
	
	return;
}
