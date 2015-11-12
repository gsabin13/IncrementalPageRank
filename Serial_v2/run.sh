#usage - run.sh 
bin=incr_pagerank
#cant.mtx  consph.mtx  mc2depi.mtx filter3D.mtx mac_econ_fwd500.mtx  mc2depi.mtx  pdb1HYS.mtx  pwtk.mtx rma10.mtx  scircuit.mtx  webbase-1M.mtx shipsec1.mtx 2cubes_sphere.mtx cage12.mtx hood.mtx m133-b3.mtx majorbasis.mtx mario002.mtx mono_500Hz.mtx offshore.mtx poisson3Da.mtx qcd5_4.mtx email-Enron.mtx facebook_combined.mtx WIPHI_graph.mtx cit-HepPh.mtx DIP_unweighted_graph.mtx  loc-gowalla_edges.mtx
# issues with patents_main.mtx cop20k_A com-amazon.ungraph.mtx roadNet-CA.mtx  twitter_combined.mtx web-NotreDame.mtx BioGRID_unweighted_graph.mtx com-dblp.ungraph.mtx  com-youtube.ungraph.mtx  web-BerkStan.mtx web-Google.mtx
for input in  pdb1HYS.mtx ;do
	mkdir -p $input
	echo $input
	
	#cp *.cu *.cuh $dir
    ./$bin /home/niuq/data/matrix_input/$input $input
done
