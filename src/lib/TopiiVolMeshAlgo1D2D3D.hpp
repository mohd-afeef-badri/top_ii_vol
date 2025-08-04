/*****************************************************************************

             This file is a part of top-ii-vol meshing tools.

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@hotmail.com
     Date     : 2019‑03‑29

     -------------------------------------------------------------------

     top-ii-vol  provides  sequential  and  parallel tools for  creating
     volumic tetrahedral meshes from a topology defined by a point cloud.
     top-ii-vol  is  distributed in the hope that it will be useful, but
     WITHOUT  ANY  WARRANTY; or  without  even  the  implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*******************************************************************************/

//-----------------------------------------------------------------------------------//
//---- Global Variables -----
//-----------------------------------------------------------------------------------//

double xx	  ;
double yy	  ;
double zz	  ;
double delz	  ;
double zznew      ;

int zglobal = pntz;
int localZmove    ;
int layerZ        ;

int IJK	          ;
int Ip1JK	  ;
int IJp1K	  ;
int IJKp1	  ;
int Ip1JKp1	  ;
int IJp1Kp1	  ;
int Ip1Jp1K	  ;
int Ip1Jp1Kp1     ;

//-----------------------------------------------------------------------------------//
//---- Setting triangle label colors -----
//-----------------------------------------------------------------------------------//

int lab_x_min = 1 ;
int lab_y_min = 2 ;
int lab_z_min = 6 ;
int lab_x_max = 4 ;
int lab_y_max = 5 ;
int lab_z_max = 3 ;

lab_y_min = 99099;
lab_y_max = 99099;
        
for(int i=0; i <NpZ; i++)
{
    for(int j=0; j <NpX; j++)
    {
        if(mpirank==i*NpX*NpY+j)
            lab_y_min = 2;
        if(mpirank==( (mpisize-1)-i*NpX*NpY  -j ))
            lab_y_max = 5;        
    }
}
 
lab_x_min = 99099;
lab_x_max = 99099;
        
for(int i=0; i <NpZ; i++)
{
    for(int j=0; j <NpY; j++)
    {
        if(mpirank==i*NpX*NpY + j*NpX )
            lab_x_min = 1;
        if(mpirank==( (mpisize-1)  -i*NpX*NpY  -j*NpX ))
            lab_x_max = 4;        
    }
}
                
lab_z_min = 99099;
lab_z_max = 99099;
        
for(int j=0; j <NpX*NpY; j++)
{
    if(mpirank==j)
        lab_z_max = 3;
    if(mpirank==( (mpisize-1)  -j))
        lab_z_min = 6;
}         
       
//-----------------------------------------------------------------------------------//
//---- Calculating Parameters -----
//-----------------------------------------------------------------------------------//

FILE* rf;
FILE* w2f ;

// Determine output format based on file extension
string outputfile_str = *outputfile+"_"+std::to_string(mpirank);
string meshtype_str = *meshtype;
if (meshtype_str == "msh" || meshtype_str == "msh2" || meshtype_str == "msh41" || meshtype_str == "msh4") {
    w2f = std::fopen((outputfile_str+".msh").c_str(), "w");
} else if (meshtype_str == "mesh") {
    w2f = std::fopen((outputfile_str+".mesh").c_str(), "w");
} else {
    std::cout << " top-ii-vol only accepts msh or mesh as mesh format " << std::endl;
}

rf = std::fopen((*inputfile+"_"+std::to_string(mpirank)+".info").c_str(), "r");
    
fscanf(rf,"%d",&pnty);
fscanf(rf,"%d",&pntx);
fscanf(rf,"%d",&pntz);
fscanf(rf,"%d",&localZmove);
fscanf(rf,"%d",&layerZ);
std::fclose(rf);

rf = std::fopen((*inputfile+"_"+std::to_string(mpirank)+".xyz").c_str(), "r");

// --------------------------------------- BUG --------------------------------------//
//                                                                                   //
//if(mpirank==(mpisize-1))pnty=pnty-1; // Last mpirank does not have comman partion  //
//                                                                                   //
// --------------------------------------- BUG --------------------------------------//

//-----------------------------------------------------------------------------------//
//---- Calculating Parameters -----
//-----------------------------------------------------------------------------------//

int PxM1 = pntx-1                                       ;
int PyM1 = pnty-1                                       ;
int PzM1 = pntz-1                                       ;
int PxPz = pntx*pntz					;
int PxPy = pntx*pnty					;

int NPnt = pntx * pnty * pntz				;
int NTri = 4*( PzM1*PxM1 + PyM1*PzM1 + PxM1*PyM1 )      ;
int NTet = PxM1 * PyM1 * PzM1 * 6			;


/*
//=============================================================================
// ------------- Commandline  output ------------------
//=============================================================================



    MPI_Barrier(MPI_COMM_WORLD);

    cout << "                                                               \n"
         << "           MPI Processs # "<<mpirank <<"                       \n"
         << "                                                               \n"
         << "  # points in the .xyz ------------- "<< PxPy << "             \n"
         << "  # points in X -------------------- "<< pntx << "             \n"
         << "  # points in Y -------------------- "<< pnty << "             \n"
         << "                                                               \n"
         << " *============================================================*\n";
*/

//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//---- Writing mesh in Medit's .mesh format -----
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//

if(meshtype_str == "mesh")
{
    //-----------------------------------------------------------------------------------//
    //---- Header for mesh -----
    //-----------------------------------------------------------------------------------//

    std::fprintf(w2f, "MeshVersionFormatted 1\n\nDimension 3\n\n");              

    //-----------------------------------------------------------------------------------//
    //---- Generating points -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout   << "Generating points ........";
    std::fprintf(w2f, "Vertices\n%d\n",NPnt);       
    for(int i=0; i<PxPy; i++)
        {
            fscanf(rf,"%lf",&xx);
            fscanf(rf,"%lf",&yy);
            fscanf(rf,"%lf",&zz);  
                  
            std::fprintf(w2f, "%lf\t%lf\t%lf 0\n"
                            ,  xx , yy ,zz
                        );

            zznew=zz;
            delz= (zmax-zz)/(((zglobal-1)-(localZmove-1*layerZ)));          
            for(int j=0; j<PzM1; j++)
                {
                    zznew  = zznew + delz;
                    std::fprintf(w2f, "%lf\t%lf\t%lf 0\n"
                                    ,  xx , yy ,zznew
                                 );                 
                }
        }

    std::fprintf(w2f, "\n");
    cout   << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 +=  t_phase;
    *time_log = string( *time_log+"\tPoint generation         : "
                            +std::to_string(t_phase)+" s\n"           );

    //-----------------------------------------------------------------------------------//
    //---- Generating Tetrahedra -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout   << "Generating Tetrahedra ...."			       ;
    std::fprintf(w2f, "Tetrahedra\n%d\n",NTet);       

    for(int j=0; j<PyM1;  j++)
        {
            for(int i=0; i<PxM1;  i++)
                {
                    for(int k=1; k<=PzM1; k++)
                        {

                            IJK	        =	i*pntz  + j*PxPz + k  ;
                            Ip1JK	    =	IJK 	+ PxPz		  ;
                            IJp1K	    =	IJK 	+ pntz		  ;
                            Ip1Jp1K     =	IJK 	+ PxPz + pntz ;
                            IJKp1       =	IJK 	+ 1			  ;
                            Ip1JKp1     =	Ip1JK 	+ 1			  ;
                            IJp1Kp1     =	IJp1K   + 1			  ;
                            Ip1Jp1Kp1   =	Ip1Jp1K + 1			  ;


                            std::fprintf(w2f, "%d\t%d\t%d\t%d 0\n"
                                              "%d\t%d\t%d\t%d 0\n"
                                              "%d\t%d\t%d\t%d 0\n"
                                              "%d\t%d\t%d\t%d 0\n"
                                              "%d\t%d\t%d\t%d 0\n"
                                              "%d\t%d\t%d\t%d 0\n"
                                            , IJK, IJKp1, IJp1K, Ip1Jp1K
                                            , IJKp1, IJK, Ip1JK, Ip1Jp1K
                                            , Ip1JKp1, IJKp1, Ip1JK, Ip1Jp1K
                                            , IJKp1, Ip1JKp1, Ip1Jp1Kp1, Ip1Jp1K
                                            , IJp1Kp1, IJKp1, Ip1Jp1Kp1, Ip1Jp1K
                                            , IJKp1, IJp1Kp1, IJp1K, Ip1Jp1K
                                         );                              
                        }
                }
        }

    std::fprintf(w2f, "\n");
    cout   << " finished for MPI rank : " << mpirank << "\n";


    t_phase = MPI_Wtime() - t_phase;
    t1 +=  t_phase;
    *time_log = string( *time_log+"\tTetrahedra generation    : "
                            +std::to_string(t_phase)+" s\n"           );

    //-----------------------------------------------------------------------------------//
    //---- Generating Triangles -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout   << "Generating Triangles ....."			       ;
    std::fprintf(w2f, "Triangles\n%d\n",NTri);       

    //----X-MIN-PLANE---//

    for(int i=0; i<PyM1;  i++)
        {
            for(int j=0; j<PzM1;  j++)
                {

                    IJK	  =	i*PxPz + j + 1	;
                    IJKp1	  =	IJK + 1         ;
                    Ip1JK	  =	IJK + PxPz      ;
                    Ip1JKp1   =	Ip1JK + 1       ;

                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                    , IJKp1, IJK, Ip1JK, lab_x_min
                                    , Ip1JKp1, IJKp1, Ip1JK, lab_x_min
                                );
                }
        }


    //----Y-MIN-PLANE----//

    for(int i=0; i<PxM1;  i++)
        {
            for(int j=0; j<PzM1;  j++)
                {

                    IJK	      =	i*pntz + j + 1	;
                    IJKp1	  =	IJK + 1		    ;
                    IJp1K	  =	IJK + pntz	    ;
                    IJp1Kp1   =	IJp1K + 1	    ;


                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                    , IJK, IJKp1, IJp1K, lab_y_min
                                    , IJKp1, IJp1Kp1, IJp1K, lab_y_min
                                );
                }
        }


    //----Z-MAX-PLANE----//

    for(int i=0; i<PyM1;  i++)
        {
            for(int j=0; j<PxM1;  j++)
                {

                    IJK	      =	i*PxPz + j*pntz +1	    ;
                    Ip1JK	  =	IJK + PxPz         		;
                    IJp1K	  =	IJK + pntz		    	;
                    Ip1Jp1K   =	Ip1JK + pntz			;

                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                    , IJK, IJp1K, Ip1Jp1K, lab_z_max
                                    , Ip1JK, IJK ,Ip1Jp1K, lab_z_max
                                );
                }
        }

    //----X-MAX-PLANE----//

    for(int i=0; i<PyM1;  i++)
        {
            for(int j=0; j<PzM1;  j++)
                {

                    IJK	      =	i*PxPz + j+1 + PxM1*(pntz)	;
                    IJKp1	  =	IJK + 1					        ;
                    Ip1JK	  =	IJK + PxPz			            ;
                    Ip1JKp1   =	Ip1JK + 1				        ;


                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                     , IJK, IJKp1, Ip1JK, lab_x_max
                                     , IJKp1, Ip1JKp1, Ip1JK, lab_x_max
                                  );
                }
        }


    //----Y-MAX-PLANE----//

    for(int i=0; i<PxM1;  i++)
        {
            for(int j=0; j<PzM1;  j++)
                {

                    IJK	      =	i*pntz + j+1 + PxPz*PyM1	;
                    IJKp1	  =	IJK + 1					    ;
                    IJp1K	  =	IJK + pntz				    ;
                    IJp1Kp1   =	IJp1K + 1				    ;

                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                    , IJKp1,IJK,IJp1K, lab_y_max
                                    , IJp1Kp1,IJKp1,IJp1K, lab_y_max
                                 );
                }
        }

    //----Z-MIN-PLANE----//

    for(int i=0; i<PyM1;  i++)
        {
            for(int j=0; j<PxM1;  j++)
                {

                    IJK	      =	i*PxPz + j*pntz + 1 + PzM1	;
                    Ip1JK	  =	IJK + PxPz			        ;
                    IJp1K	  =	IJK + pntz				    ;
                    Ip1Jp1K   =	Ip1JK + pntz				;

                    std::fprintf(w2f, "%d\t%d\t%d\t%d\n"
                                      "%d\t%d\t%d\t%d\n"
                                    , IJp1K, IJK, Ip1Jp1K, lab_z_min
                                    , IJK, Ip1JK, Ip1Jp1K, lab_z_min
                                 );
                }
        }

    cout   << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 +=  t_phase;
    *time_log = string( *time_log+"\tTriangles generation     : "
                            +std::to_string(t_phase)+" s\n"           );

    //-----------------------------------------------------------------------------------//
    //---- Finishing footer -----
    //-----------------------------------------------------------------------------------//
           
    std::fprintf(w2f, "\nEnd\n");
}

//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//---- Writing mesh in Gmsh's .msh 4.1 format -----
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//

else if (meshtype_str == "msh4" || meshtype_str == "msh41")
{
    cout << "Meshing the topology in Gmsh's *.msh 4.1 format for MPI rank : " << mpirank << "\n";

    //--------------------------------------//
    //---- Header for msh4.1 -----
    //--------------------------------------//

    std::fprintf(w2f, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");

    // Physical Names section
    std::fprintf(w2f, "$PhysicalNames\n");
    std::fprintf(w2f, "7\n");              // number of physical groups

    std::fprintf(w2f, "2 %d \"%s\"\n", lab_x_min, (lab_x_min != 99099 ? "left" : "interface"));   // X-min surface
    std::fprintf(w2f, "2 %d \"%s\"\n", lab_y_min, (lab_y_min != 99099 ? "front" : "interface"));  // Y-min surface
    std::fprintf(w2f, "2 %d \"%s\"\n", lab_z_max, (lab_z_max != 99099 ? "top" : "interface"));    // Z-max surface
    std::fprintf(w2f, "2 %d \"%s\"\n", lab_x_max, (lab_x_max != 99099 ? "right" : "interface"));  // X-max surface
    std::fprintf(w2f, "2 %d \"%s\"\n", lab_y_max, (lab_y_max != 99099 ? "back" : "interface"));   // Y-max surface
    std::fprintf(w2f, "2 %d \"%s\"\n", lab_z_min, (lab_z_min != 99099 ? "bottom" : "interface")); // Z-min surface

    std::fprintf(w2f, "3 7 \"volume\"\n");             // Volume (physical tag 7)
    std::fprintf(w2f, "$EndPhysicalNames\n");

    // Entities section (required in MSH 4.1)
    std::fprintf(w2f, "$Entities\n");
    // Format: numPointEntities numCurveEntities numSurfaceEntities numVolumeEntities
    std::fprintf(w2f, "0 0 6 1\n");

    // Define 6 surfaces (boundary faces) with proper physical tags
    // Format: surfaceTag minX minY minZ maxX maxY maxZ numBoundingCurves physicalTag
    std::fprintf(w2f, "1 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_x_min); // X-min surface
    std::fprintf(w2f, "2 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_y_min); // Y-min surface
    std::fprintf(w2f, "3 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_z_max); // Z-max surface
    std::fprintf(w2f, "4 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_x_max); // X-max surface
    std::fprintf(w2f, "5 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_y_max); // Y-max surface
    std::fprintf(w2f, "6 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 %d 0\n", lab_z_min); // Z-min surface

    // Define 1 volume
    // Format: volumeTag minX minY minZ maxX maxY maxZ numBoundingSurfaces physicalTag
    std::fprintf(w2f, "1 -1e99 -1e99 -1e99 1e99 1e99 1e99 1 7 0\n");

    std::fprintf(w2f, "$EndEntities\n");

    // Nodes section
    std::fprintf(w2f, "$Nodes\n");
    std::fprintf(w2f, "1 %d 1 %d\n", NPnt, NPnt); // numEntityBlocks numNodes minNodeTag maxNodeTag
    std::fprintf(w2f, "3 1 0 %d\n", NPnt);        // entityDim entityTag parametric numNodesInBlock

    //-----------------------------------------------------------------------------------//
    //---- Generating points -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout << "Generating points ........";

    // Node tags
    for (int i = 1; i <= NPnt; i++)
    {
        std::fprintf(w2f, "%d\n", i);
    }

    int counter1 = 1;

    for (int i = 0; i < PxPy; i++)
    {
        fscanf(rf, "%lf", &xx);
        fscanf(rf, "%lf", &yy);
        fscanf(rf, "%lf", &zz);

        std::fprintf(w2f, "%lf %lf %lf\n", xx, yy, zz);

        counter1++;
        zznew = zz;
        delz = (zmax - zz) / (((zglobal - 1) - (localZmove - 1 * layerZ)));
        for (int j = 0; j < PzM1; j++)
        {
            zznew = zznew + delz;
            std::fprintf(w2f, "%lf %lf %lf\n", xx, yy, zznew);
            counter1++;
        }
    }

    std::fprintf(w2f, "$EndNodes\n");
    cout << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 += t_phase;
    *time_log = string(*time_log + "\tPoint generation         : " +
                       std::to_string(t_phase) + " s\n");

    //-----------------------------------------------------------------------------------//
    //---- Generating Elements (Triangles and Tetrahedra) -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout << "Generating Elements ....";
    std::fprintf(w2f, "$Elements\n");

    // Calculate number of element blocks (6 surface blocks + 1 volume block)
    int numElementBlocks = 7;
    std::fprintf(w2f, "%d %d 1 %d\n", numElementBlocks, NTet + NTri, NTet + NTri);

    counter1 = 1;

    //---------------------------------X-MIN-PLANE (Surface 1)----------------------------------//
    int numTriXMin = 2 * PyM1 * PzM1;
    std::fprintf(w2f, "2 1 2 %d\n", numTriXMin); // entityDim entityTag elementType numElementsInBlock

    for (int i = 0; i < PyM1; i++)
    {
        for (int j = 0; j < PzM1; j++)
        {
            IJK = i * PxPz + j + 1;
            IJKp1 = IJK + 1;
            Ip1JK = IJK + PxPz;
            Ip1JKp1 = Ip1JK + 1;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJKp1, IJK, Ip1JK);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, Ip1JKp1, IJKp1, Ip1JK);
            counter1 += 2;
        }
    }

    //---------------------------------Y-MIN-PLANE (Surface 2)----------------------------------//
    int numTriYMin = 2 * PxM1 * PzM1;
    std::fprintf(w2f, "2 2 2 %d\n", numTriYMin);

    for (int i = 0; i < PxM1; i++)
    {
        for (int j = 0; j < PzM1; j++)
        {
            IJK = i * pntz + j + 1;
            IJKp1 = IJK + 1;
            IJp1K = IJK + pntz;
            IJp1Kp1 = IJp1K + 1;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJK, IJKp1, IJp1K);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, IJKp1, IJp1Kp1, IJp1K);
            counter1 += 2;
        }
    }

    //---------------------------------Z-MAX-PLANE (Surface 3)----------------------------------//
    int numTriZMax = 2 * PyM1 * PxM1;
    std::fprintf(w2f, "2 3 2 %d\n", numTriZMax);

    for (int i = 0; i < PyM1; i++)
    {
        for (int j = 0; j < PxM1; j++)
        {
            IJK = i * PxPz + j * pntz + 1;
            Ip1JK = IJK + PxPz;
            IJp1K = IJK + pntz;
            Ip1Jp1K = Ip1JK + pntz;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJK, IJp1K, Ip1Jp1K);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, Ip1JK, IJK, Ip1Jp1K);
            counter1 += 2;
        }
    }

    //---------------------------------X-MAX-PLANE (Surface 4)----------------------------------//
    int numTriXMax = 2 * PyM1 * PzM1;
    std::fprintf(w2f, "2 4 2 %d\n", numTriXMax);

    for (int i = 0; i < PyM1; i++)
    {
        for (int j = 0; j < PzM1; j++)
        {
            IJK = i * PxPz + j + 1 + PxM1 * pntz;
            IJKp1 = IJK + 1;
            Ip1JK = IJK + PxPz;
            Ip1JKp1 = Ip1JK + 1;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJK, IJKp1, Ip1JK);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, IJKp1, Ip1JKp1, Ip1JK);
            counter1 += 2;
        }
    }

    //---------------------------------Y-MAX-PLANE (Surface 5)----------------------------------//
    int numTriYMax = 2 * PxM1 * PzM1;
    std::fprintf(w2f, "2 5 2 %d\n", numTriYMax);

    for (int i = 0; i < PxM1; i++)
    {
        for (int j = 0; j < PzM1; j++)
        {
            IJK = i * pntz + j + 1 + PxPz * PyM1;
            IJKp1 = IJK + 1;
            IJp1K = IJK + pntz;
            IJp1Kp1 = IJp1K + 1;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJKp1, IJK, IJp1K);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, IJp1Kp1, IJKp1, IJp1K);
            counter1 += 2;
        }
    }

    //---------------------------------Z-MIN-PLANE (Surface 6)----------------------------------//
    int numTriZMin = 2 * PyM1 * PxM1;
    std::fprintf(w2f, "2 6 2 %d\n", numTriZMin);

    for (int i = 0; i < PyM1; i++)
    {
        for (int j = 0; j < PxM1; j++)
        {
            IJK = i * PxPz + j * pntz + 1 + PzM1;
            Ip1JK = IJK + PxPz;
            IJp1K = IJK + pntz;
            Ip1Jp1K = Ip1JK + pntz;

            std::fprintf(w2f, "%d %d %d %d\n", counter1, IJp1K, IJK, Ip1Jp1K);
            std::fprintf(w2f, "%d %d %d %d\n", counter1 + 1, IJK, Ip1JK, Ip1Jp1K);
            counter1 += 2;
        }
    }

    //-----------------------------------------------------------------------------------//
    //---- Generating Tetrahedra (Volume 1) -----
    //-----------------------------------------------------------------------------------//

    cout << " finished triangles for MPI rank : " << mpirank << "\n";
    cout << "Generating Tetrahedrals ....";

    // Volume elements
    std::fprintf(w2f, "3 1 4 %d\n", NTet); // entityDim entityTag elementType(4=tetrahedra) numElementsInBlock

    for (int j = 0; j < PyM1; j++)
    {
        for (int i = 0; i < PxM1; i++)
        {
            for (int k = 1; k <= PzM1; k++)
            {
                IJK = i * pntz + j * PxPz + k;
                Ip1JK = IJK + PxPz;
                IJp1K = IJK + pntz;
                Ip1Jp1K = IJK + PxPz + pntz;
                IJKp1 = IJK + 1;
                Ip1JKp1 = Ip1JK + 1;
                IJp1Kp1 = IJp1K + 1;
                Ip1Jp1Kp1 = Ip1Jp1K + 1;

                std::fprintf(w2f, "%d %d %d %d %d\n", counter1, IJK, IJKp1, IJp1K, Ip1Jp1K);
                std::fprintf(w2f, "%d %d %d %d %d\n", counter1 + 1, IJKp1, IJK, Ip1JK, Ip1Jp1K);
                std::fprintf(w2f, "%d %d %d %d %d\n", counter1 + 2, Ip1JKp1, IJKp1, Ip1JK, Ip1Jp1K);
                std::fprintf(w2f, "%d %d %d %d %d\n", counter1 + 3, IJKp1, Ip1JKp1, Ip1Jp1Kp1, Ip1Jp1K);
                std::fprintf(w2f, "%d %d %d %d %d\n", counter1 + 4, IJp1Kp1, IJKp1, Ip1Jp1Kp1, Ip1Jp1K);
                std::fprintf(w2f, "%d %d %d %d %d\n", counter1 + 5, IJKp1, IJp1Kp1, IJp1K, Ip1Jp1K);

                counter1 += 6;
            }
        }
    }

    std::fprintf(w2f, "$EndElements\n");
    cout << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 += t_phase;
    *time_log = string(*time_log + "\tElements generation      : " +
                       std::to_string(t_phase) + " s\n");
}

//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//---- Writing mesh in Gmsh's .msh format -----
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//

else if (meshtype_str == "msh")
{
    cout << "Meshing the topology in Gmsh's *.msh format for MPI rank : " << mpirank << "\n";

    //--------------------------------------//
    //---- Header for msh2 -----
    //--------------------------------------//

    std::fprintf(w2f, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n",NPnt);

    int counter1 = 1;

    //-----------------------------------------------------------------------------------//
    //---- Generating points -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout << "Generating points ........";

    for(int i=0; i<PxPy; i++)
    {
        fscanf(rf,"%lf",&xx);
        fscanf(rf,"%lf",&yy);
        fscanf(rf,"%lf",&zz); 

        std::fprintf(w2f, "%d %lf %lf %lf\n"
                        ,  counter1,xx,yy,zz
                    );

        counter1++ ;
        zznew = zz ;
        delz= (zmax-zz)/(((zglobal-1)-(localZmove-1*layerZ)));
        for(int j=0; j<PzM1; j++)
        {
            zznew  = zznew + delz ;
            std::fprintf(w2f, "%d %lf %lf %lf\n"
                            ,  counter1,xx,yy,zznew
                        );
            counter1++;
        }
    }

    std::fprintf(w2f, "$EndNodes\n");
    cout << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 +=  t_phase;
    *time_log = string( *time_log+"\tPoint generation         : "
                            +std::to_string(t_phase)+" s\n"           );

    //-----------------------------------------------------------------------------------//
    //---- Generating Triangles and Tetrahedra -----
    //-----------------------------------------------------------------------------------//

    t_phase = MPI_Wtime();

    cout << "Generating Triangles and Tetrahedra ....";
    std::fprintf(w2f, "$Elements\n%d\n",NTet+NTri);
    
    counter1=1;

    //---------------------------------X-MIN-PLANE----------------------------------//
    for(int i=0; i<PyM1;  i++)
    {
        for(int j=0; j<PzM1;  j++)
        {
            IJK	    = i*PxPz + j + 1;
            IJKp1   = IJK + 1;
            Ip1JK   = IJK + PxPz;
            Ip1JKp1 = Ip1JK + 1;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_x_min, lab_x_min*1000+lab_x_min, IJKp1, IJK, Ip1JK
                            , counter1+1, lab_x_min, lab_x_min*1000+lab_x_min, Ip1JKp1, IJKp1, Ip1JK
                        );

            counter1=counter1+2;
        }
    }

    //---------------------------------Y-MIN-PLANE----------------------------------//
    for(int i=0; i<PxM1;  i++)
    {
        for(int j=0; j<PzM1;  j++)
        {
            IJK	    = i*pntz + j + 1;
            IJKp1   = IJK + 1;
            IJp1K   = IJK + pntz;
            IJp1Kp1 = IJp1K + 1;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_y_min, lab_y_min*1000+lab_y_min, IJK, IJKp1, IJp1K
                            , counter1+1, lab_y_min, lab_y_min*1000+lab_y_min, IJKp1, IJp1Kp1, IJp1K
                        );

            counter1=counter1+2;
        }
    }

    //---------------------------------Z-MAX-PLANE----------------------------------//
    for(int i=0; i<PyM1;  i++)
    {
        for(int j=0; j<PxM1;  j++)
        {
            IJK	    = i*PxPz + j*pntz + 1;
            Ip1JK   = IJK + PxPz;
            IJp1K   = IJK + pntz;
            Ip1Jp1K = Ip1JK + pntz;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_z_max, lab_z_max*1000+lab_z_max, IJK, IJp1K, Ip1Jp1K
                            , counter1+1, lab_z_max, lab_z_max*1000+lab_z_max, Ip1JK, IJK, Ip1Jp1K
                        );

            counter1=counter1+2;
        }
    }

    //---------------------------------X-MAX-PLANE----------------------------------//
    for(int i=0; i<PyM1;  i++)
    {
        for(int j=0; j<PzM1;  j++)
        {
            IJK	    = i*PxPz + j + 1 + PxM1*pntz;
            IJKp1   = IJK + 1;
            Ip1JK   = IJK + PxPz;
            Ip1JKp1 = Ip1JK + 1;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_x_max, lab_x_max*1000+lab_x_max, IJK, IJKp1, Ip1JK
                            , counter1+1, lab_x_max, lab_x_max*1000+lab_x_max, IJKp1, Ip1JKp1, Ip1JK
                        );
            
            counter1=counter1+2;
        }
    }

    //---------------------------------Y-MAX-PLANE----------------------------------//
    for(int i=0; i<PxM1;  i++)
    {
        for(int j=0; j<PzM1;  j++)
        {
            IJK	    = i*pntz + j + 1 +  PxPz*PyM1;
            IJKp1   = IJK + 1;
            IJp1K   = IJK + pntz;
            IJp1Kp1 = IJp1K + 1;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_y_max, lab_y_max*1000+lab_y_max, IJKp1, IJK, IJp1K
                            , counter1+1, lab_y_max, lab_y_max*1000+lab_y_max, IJp1Kp1, IJKp1, IJp1K
                        );
            
            counter1=counter1+2;
        }
    }

    //---------------------------------Z-MIN-PLANE----------------------------------//
    for(int i=0; i<PyM1;  i++)
    {
        for(int j=0; j<PxM1;  j++)
        {
            IJK     = i*PxPz + j*pntz + 1 + PzM1;
            Ip1JK   = IJK + PxPz;
            IJp1K   = IJK + pntz;
            Ip1Jp1K = Ip1JK + pntz;

            std::fprintf(w2f, "%d 2 2 %d %d %d %d %d\n"
                              "%d 2 2 %d %d %d %d %d\n"
                            , counter1, lab_z_min, lab_z_min*1000+lab_z_min, IJp1K, IJK, Ip1Jp1K
                            , counter1+1, lab_z_min, lab_z_min*1000+lab_z_min, IJK, Ip1JK, Ip1Jp1K
                        );

            counter1 = counter1+2;
        }
    }

    //-----------------------------------------------------------------------------------//
    //---- Generating Tetrahedra -----
    //-----------------------------------------------------------------------------------//

    cout << " finished triangles for MPI rank : " << mpirank << "\n";
    cout << "Generating Tetrahedrals ....";

    for(int j=0; j<PyM1;  j++)
    {
        for(int i=0; i<PxM1;  i++)
        {
            for(int k=1; k<=PzM1; k++)
            {
                IJK	      =	i*pntz + j*PxPz + k;
                Ip1JK     =	IJK + PxPz;
                IJp1K     =	IJK + pntz;
                Ip1Jp1K   =	IJK + PxPz + pntz;
                IJKp1     =	IJK + 1;
                Ip1JKp1   =	Ip1JK + 1;
                IJp1Kp1   =	IJp1K + 1;
                Ip1Jp1Kp1 =	Ip1Jp1K + 1;

                std::fprintf(w2f, "%d 4 2 7 77 %d %d %d %d\n"
                                  "%d 4 2 7 77 %d %d %d %d\n"
                                  "%d 4 2 7 77 %d %d %d %d\n"
                                  "%d 4 2 7 77 %d %d %d %d\n"
                                  "%d 4 2 7 77 %d %d %d %d\n"
                                  "%d 4 2 7 77 %d %d %d %d\n"
                               , counter1  , IJK, IJKp1, IJp1K,Ip1Jp1K
                               , counter1+1, IJKp1, IJK, Ip1JK, Ip1Jp1K
                               , counter1+2, Ip1JKp1, IJKp1, Ip1JK, Ip1Jp1K
                               , counter1+3, IJKp1, Ip1JKp1, Ip1Jp1Kp1, Ip1Jp1K
                               , counter1+4, IJp1Kp1, IJKp1, Ip1Jp1Kp1, Ip1Jp1K
                               , counter1+5, IJKp1, IJp1Kp1, IJp1K, Ip1Jp1K
                            );

                counter1=counter1+6;
            }
        }
    }

    std::fprintf(w2f, "$EndElements\n");
    cout  << " finished for MPI rank : " << mpirank << "\n";

    t_phase = MPI_Wtime() - t_phase;
    t1 +=  t_phase;
    *time_log = string( *time_log+"\tTriangles & Tetrahedra generation : "
                            +std::to_string(t_phase)+" s\n"           );
}

std::fclose(w2f);
std::fclose(rf);