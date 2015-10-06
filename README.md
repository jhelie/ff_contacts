# ff_contacts
Python script to calculate contacts between flip-flopping lipids and other beads in the system

###DESCRIPTION


This script produces contact statistics between flipflopping lipids and transmembrane
protein clusters.

A file listing the flip-flopping lipids must be supplied with the --flipflops option.
Each line of this file should follow the format (time in ns):

 -> 'resname,resid,starting_leaflet,z_bead,t_start,t_end'

where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4,150,500'.
The 'z_bead' particle is used to track the position of the lipid. Such a file can be
obtain by running successively the ff_detect and ff_times scripts.

This script produced 2 types of ouptus::
 -> contact statistics between flip-flopping lipids and TM clusters
 -> with --profile: distribution profile of those contacts along the local bilayer normal

#####Detection of transmembrane protein clusters
-------------------------------------------
Two clustering algorithms can be used to identify protein clusters.
->Connectivity based (relies on networkX module):
  A protein is considered in a cluster if it is within a distance less than --nx_cutoff
  from another protein. This means that a single protein can act as a connector between
  two otherwise disconnected protein clusters.
  This algorithm can be ran using either the minimum distante between proteins (default, 
  --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
  The 'min' option scales as the square of the number of proteins and can thus be very
  slow for large systems.

->Density based (relies on the sklearn module and its implementation of DBSCAN):
  A protein is considered in a cluster if is surrounded by at least --db_neighbours other
  proteins within a radius of --db_radius.
  This density based approach is usually less suited to the detection of protein
  clusters but as a general rule the more compact the clusters, the smaller --db_radius
  the higher --db_neighbours can be - for details on this algorithm see its online
  documentation.
  This algorithm is selected by setting the --algorithm option to 'density'.

The identified protein clusters are considered to be transmembrane only if the closest
lipid headgroup neighbours to the cluster particles are all within the same leaflet.
In addition to the sizes identified, size groups can be defined - see note 4(c).


###REQUIREMENTS

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - networkX (if --algorithm set to 'min' or 'cog')
 - sklearn (if --algorithm set to 'density')


###NOTES

1. The types of contacts to take into account must be specified via --contacts. All
   unecessary particles (eg water particles and lipids tails) should be removed from
   the trajectory for the script to run as fast as possible.

2. Flip-flopping lipids and their start/end times can be identified by successively
   running the ff_detect and ff_times scripts.

3. Identification of the bilayer leaflets can be controlled via two options.
   (a) beads
    By default, the particles taken into account to define leaflet depend on the
    forcefield (which can be set via the --forcefield option) and are as follows:
    -> Martini: 'name PO4 or name PO3 or name B1A'
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
     -> '--leaflet 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflet large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the 1st frame of the xtc
    file supplied in order to get a meaningful outcome. 

	NOTE: By default the gro file is only used as a topology file and the 1st frame of the
	xtc is used to identify leaflets. If you wish to use the gro file instead, for instance
	in the case that the 1st frame of the xtc is not flat, you need to specify the --use_gro
	flag: be warned that this might take a few minutes longer on large systems.

4. (a) Proteins are detected automatically but you can specify an input file to define your
   own selection with the --proteins option.
   In this case the supplied file should contain on each line a protein selection string
   that can be passed as the argument of the MDAnalysis selectAtoms() routine - for 
   instance 'bynum 1:344'.

   (b) The size detected for a TM cluster can vary slightly between frames (due to
   peripheral or 'loosely' attached peptide moving in and out the cutoff). The detected
   size of a cluster is simply recorded and the data is just binned into each size or size
   group without clever analysis (it would be futile to attempt to define what is 'the'
   size of a dynamic and transient cluster).
   The number of sampling events for which each cluster size/group has been sampled is
   presented in the outputs so that you can judge how representative each size/group
   actually is.

   (c) Cluter sizes groups can ebe specified via --groups to bin conctact statistics.
   A weighted average of individual cluster size data, using the number of events
   sampled, is used to calculate groups statistics. Even when size groups are specified,
   individual data for each size detected is also  generated. The sizes groups are defined
   by supplying a file with the --groups option, whose lines all follow the format:
    -> 'lower_size,upper_size, colour'

   Size groups definition should follow the following rules:
    -to specify an open ended group use 'max', e.g. '3,max'
    -groups should be ordered by increasing size and their boundaries should not overlap
    -boundaries are inclusive so you can specify unique size groups with 'size,size'
    -any cluster size not falling within the specified size groups will be labeled as 'other'

5. There are 3 possible options to determine the local normal to the bilayer. These are
   controlled with the flags --normal and --normal_d:
   (a) 'z': the bilayer is assumed flat in the x,y plane and the z axis is taken to be the
    normal. Best for systems without significant curvature and local deformations. In this
    case the --normal_d flag is ignored.

   (b) 'cog': in this case neighbourhing particles to current cluster of interest are
    identified in the lower and upper leaflet. The local normal is then considered to be the
    vector going from the cog of the lower ones to the cog of the upper ones. In each leaflet,
    neighbouring particles are the particles selected by --beads which are within --normal_d
    Angstrom of the cog of the protein cluster of interest.

   (c) 'svd': in this case neighbourhing particles to current cluster of interest are
    identified in the lower and upper leaflet as in (b) above. The normal of the best fitting
    plane to these particles is obtained by performing a singular value decomposition of their
    coordinate matrix.


###USAGE
```	
Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		1	: process every t-frames

Lipids identification  
-----------------------------------------------------
--flipflops		: input file with flipflopping lipids, see note 2
--beads			: leaflet identification technique, see note 3(a)
--leaflets	optimise: leaflet identification technique, see note 3(b)
--use_gro		: use gro file instead of xtc, see note 3(b)

Protein clusters identification and contacts
-----------------------------------------------------
--proteins		: protein selection file, (optional, see note 6)
--groups		: cluster groups definition file, see note 4(c)
--algorithm	min	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 	6	: networkX cutoff distance for protein-protein contact (Angstrom)
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	

 Contacts profile options
-----------------------------------------------------
--profile 		: calculate local bilayer normal and store position of each contacts
--normal	svd	: local normal to bilayer ('z', 'cog' or 'svd'), see note 5
--normal_d	50	: distance of points to take into account for local normal, see note 5
--pl_cutoff 	6	: cutoff distance for protein-lipid contact (Angstrom)
--slices_dist	40 	: max distance from center of bilayer (Angstrom)
--slices_thick	0.5 	: thickness of the profile slices (Angstrom)
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
```
