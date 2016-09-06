//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_2D.h
//OBJECTIVE:	Create a 2D nanotube netwok
//AUTHOR:		Fei Han; Maloth Thirupathi
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef APPNETWORK2D_H
#define APPNETWORK2D_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "GenNetwork_2D.h"
#include "Cutoff_Wins.h"
#include "Hoshen_Kopelman.h"
//#include "Direct_Electrifying.h"
//#include "Backbone_Network.h"
#include "Background_vectors.h"
#include "Contact_grid.h"
#include "Percolation.h"
//#include "Clusters_fractions.h"
#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
class App_Network
{
public:
    //Data Member
    //vector<Point_2D> cnps;			//Define 3D point verctor of nanotuber points
    //vector<double> cnts_radius;		//Define the radius of every nanotube in the network
    
    //Constructor
    App_Network(){};
    
    //Member Functions
    int Create_conductive_network_2D(Input *Init)const;
    //int Export_tecplot_files(const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &isolated, vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_indices)const;
    //int Convert_index_to_structure(const vector<long int> &indices, vector<vector<long int> > &structure)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================

