/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include "hashtables_bfields.h"
#include "stats.h"
#include "hashmap.h"
#include "tree.h"
#include "tree_utils.h"

/* Returns a table of all node ids of the tree, with 1 if they are taxon on the side of the edge, 0 if not (or internal) */
int fill_all_taxa_ids(Node *node, Node *prev, int *output){
  int i;
  int nbt = 0;
  if(node->nneigh == 1){
    output[node->id] = 1;
    return 1;
  } else{
    for(i=0; i < node->nneigh; i++){
      if(node->neigh[i] != prev)
	nbt+=fill_all_taxa_ids(node->neigh[i], node, output);
    }
    return(nbt);
  }
}


/* Returns a table of a sample of the node ids that are taxon and on the orientation given */
int * sample_taxa(Tree *t, int n_to_sample, Node *node, Node *prev){
  int * allnodes = (int*) calloc(t->nb_nodes, sizeof(int));
  int nbtax = fill_all_taxa_ids(node,prev,allnodes);
  int * alltax = (int*) calloc(nbtax, sizeof(int));
  int *output;
  int cur = 0;
  int i;
  for(i=0; i < t->nb_nodes; i++){
    if(allnodes[i]){
      alltax[cur] = i;
      cur++;
    }
  }
  output = sample(alltax, nbtax, n_to_sample, 0);
  free(allnodes);
  free(alltax);

  return(output);
}

/**
   This method swaps 2 edges connected to the given edge.
   If e is terminal, does nothing
   Before
     a       d 
      \  e  /  
   left.---.right
      /     \	 
     b       c 

   After
    c       d        a       b
     \     /  	      \     /  
      .---.     or     .---.    
     /     \	      /     \	 
    b       a 	     d       c 

    Randomly one of the two options
    returns the min topo_depth of the 2 swaped branches
 */
int swap_branches(Tree *t, Edge *e){
  Node *left = e->left;
  Node *right = e->right;
  if(left->nneigh == 1 || right->nneigh==1){
    return(0);
  }

  int dir_left_to_right = dir_a_to_b(left, right);
  int dir_right_to_left = dir_a_to_b(right, left);

  /* The two edges are chosen randomly in the left and in the right of e */
  int picked_left_index = rand_to(left->nneigh-1)+1;
  int picked_right_index = rand_to(right->nneigh-1)+1;

  picked_left_index = (dir_left_to_right+picked_left_index)%left->nneigh;
  picked_right_index= (dir_right_to_left+picked_right_index)%right->nneigh;

  Edge *picked_left_branch =  left->br[picked_left_index];
  Edge *picked_right_branch = right->br[picked_right_index];

  int i;
  /* fprintf(stderr,"left  branch %d | Topo= %d\n",picked_left_branch->id,picked_left_branch->topo_depth); */
  /* fprintf(stderr,"right branch %d | Topo= %d\n",picked_right_branch->id,picked_right_branch->topo_depth); */
  /* for(i=0;i<t->nb_taxa;i++){ */
  /*   if(lookup_id(t->a_edges[picked_left_branch->id]->hashtbl[1],i)) */
  /*     fprintf(stderr," %s",t->taxa_names[i]); */
  /* } */
  /* fprintf(stderr,"\n"); */
  /* for(i=0;i<t->nb_taxa;i++){ */
  /*   if(lookup_id(t->a_edges[picked_right_branch->id]->hashtbl[1],i)) */
  /*     fprintf(stderr," %s",t->taxa_names[i]); */
  /* } */
  /* fprintf(stderr,"\n"); */


  int sum_depth = picked_left_branch->topo_depth + picked_right_branch->topo_depth;

  /**
     All the other participants to this swap (see figure)
   */
  Node *a,*c;
  int a_to_left_dir,
    left_to_a_dir;
  int c_to_right_dir,
    right_to_c_dir;

  if(picked_left_branch->right==left){
    a=picked_left_branch->left;
  }else{
    a=picked_left_branch->right;
  }
  if(picked_right_branch->right==right){
    c=picked_right_branch->left;
  }else{
    c=picked_right_branch->right;
  }
  a_to_left_dir = dir_a_to_b(a, left);
  left_to_a_dir = dir_a_to_b(left, a);
  c_to_right_dir = dir_a_to_b(c, right);
  right_to_c_dir = dir_a_to_b(right, c);

  /**
     We swap the two edges 
  */

  /* First swap the edge pointers of the edges */
  if(picked_left_branch->right==left){
    picked_left_branch->right = right;
  }else{
    picked_left_branch->left = right;
  }
  if(picked_right_branch->right==right){
    picked_right_branch->right = left;
  }else{
    picked_right_branch->left = left;
  }
  
  /*We then swap the node pointers of the nodes*/
  a->neigh[a_to_left_dir] = right;
  c->neigh[c_to_right_dir]= left;
  left->neigh[left_to_a_dir] = c;
  right->neigh[right_to_c_dir] = a;

  /* And the final edges pointers of the nodes */
  left->br[picked_left_index] = picked_right_branch;
  right->br[picked_right_index] = picked_left_branch;

  /**
     We recompute hashtables and node depths
   */
  for (i = 0; i < t->nb_edges; i++) {
    if(t->a_edges[i]->hashtbl[0] != NULL)
      free_id_hashtable(t->a_edges[i]->hashtbl[0]);
    if(t->a_edges[i]->hashtbl[1] != NULL)
      free_id_hashtable(t->a_edges[i]->hashtbl[1]); 
    t->a_edges[i]->hashtbl[0] = create_id_hash_table(t->length_hashtables);
    t->a_edges[i]->hashtbl[1] = create_id_hash_table(t->length_hashtables);
  }
 
  update_hashtables_post_alltree(t);
  update_hashtables_pre_alltree(t);
  update_node_depths_post_alltree(t);
  update_node_depths_pre_alltree(t);

  for (i = 0; i < t->nb_edges; i++) {
    free_id_hashtable(t->a_edges[i]->hashtbl[0]); 
    t->a_edges[i]->hashtbl[0] = NULL;
  }

  /* topological depths of branches */
  update_all_topo_depths_from_hashtables(t);
  return(sum_depth);
}

/**
   Here we will test the classical bootstrap with a very simple case (to test hashtables)
   
 */
int test_classical_bootstrap(){
  /*
    Tree 1          Bootstrap tree: 
     a        d     a        d   
      \  e   /       \   e  /   
       .---(.)	     (.)---.	     
      /      \	     /      \	     
     b        c     b        c  
      The node (.) is the top node of the newick file: 
      It changes the orientation of the hashtables for the edge e
      2 newick representations of the SAME tree
   */
  char *ref_tree_string   = "((a:1,b:1):1,c:1,d:1);";
  char *boot_tree_string = "(b:1,a:1,(c:1,d:1):1);";
  char** taxname_lookup_table = NULL;

  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  Tree* boot_tree = complete_parse_nh(boot_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  int i,j;
  int common_splits = 0;
  int splits_not_found = 0;
  for (i = 0; i < ref_tree->nb_edges; i++) {
    if(ref_tree->a_edges[i]->right->nneigh == 1) continue;
    /* we skip the branches leading to leaves */
    
    if(ref_tree->a_edges[i]->had_zero_length) continue;
    /* a branch with length == 0 is not to be considered as a valid bipartition */
    
    for (j = 0; j < boot_tree->nb_edges; j++) {
      if(boot_tree->a_edges[j]->had_zero_length) continue;
      /* a branch with length == 0 is not to be considered as a valid bipartition */
      if (equal_or_complement_id_hashtables(ref_tree->a_edges[i]->hashtbl[1],
					    boot_tree->a_edges[j]->hashtbl[1],
					    ref_tree->nb_taxa)) {
	//printf("result: splits ARE equal!\n");
	common_splits++;
	break;
      }
    } /* end for on j */
    if (j == boot_tree->nb_edges) splits_not_found++;
  } /* end for on i */
  free_tree(boot_tree);
  free_tree(ref_tree);
  free(taxname_lookup_table); /* which is a (char**) */

  if(common_splits != 1){
    fprintf(stderr,"Classical Bootstrap test error: Number of common splits is: %d, and should be: %d\n",common_splits,1);
    return EXIT_FAILURE;
  }
  if(splits_not_found != 0){
    fprintf(stderr,"Classical Bootstrap test error: Number of splits not found is: %d, and should be: %d\n",splits_not_found,0);
    return EXIT_FAILURE;
  }
  fprintf(stderr,"Classical Bootstrap test : OK\n");
  return(EXIT_SUCCESS);
}

void test_fill_hashtable_post_order(Node* current, Node* orig, Tree* t, id_hash_table_t *h) {
	/* we are going to update one of the two hashtables sitting on the branch between current and orig. */
	int i, n = current->nneigh;
	if(orig == NULL) return;
	int curr_to_orig = dir_a_to_b(current, orig);

	Edge* br = current->br[curr_to_orig]; /* br: current to orig; br2: any _other_ branch from current */

	for(i=1 ; i < n ; i++) {
	  test_fill_hashtable_post_order(current->neigh[(curr_to_orig+i)%n], current,t,h);
	}

	/* but if n = 1 we haven't done anything (leaf): we must put the info corresponding to the taxon into the branch */
	if (n == 1) {
	  assert(br->right == current);
	  /* add the id of the taxon to the right hashtable of the branch */
	  add_id(h,get_tax_id_from_tax_name(current->name, t->taxname_lookup_table, t->nb_taxa));
	}
} /* end update_hashtables_post_doer */




int test_swap_branches(){
  /**
      a     e     d 
       \    |    /  
        .---.---.   
       /  *   *  \	 
      b           c 

      We will swap the edges from one of the edges * 
   */
  char** taxname_lookup_table = NULL;
  char *ref_tree_string = "((a:1,b:1):1,e:1,(c:1,d:1):1);"; 
  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */

  int e_index;
  int swaped = 0;
  for(e_index=0;e_index<ref_tree->nb_edges;e_index++){
    Edge *e = ref_tree->a_edges[e_index];
    if(e->right->nneigh>1 &&
       e->left->nneigh>1 
       && !swaped){
      /*On swap la première qui vient*/
      swap_branches(ref_tree,e);
      swaped = 1;
    }
  }
  fprintf(stderr,"Swap branch Test: OK\n");
  return(EXIT_SUCCESS);
}

/**
   We test the TRANSFER Support for branches of the initial tree compared to another tree

 */
int test_transfer_1(){
  srand(time(NULL));
  char *ref_tree_string = "((a:1,b:1,c:1):1,(d:1,e:1,f:1):1,((g:1,h:1,i:1):1,(j:1,k:1,l:1):1,(m:1,n:1,o:1):1):1);"; 
  char *swap_tree_string = "((g:1.000000,h:1.000000,i:1.000000):1.000000,(m:1.000000,n:1.000000,o:1.000000):1.000000,((a:1.000000,b:1.000000,c:1.000000):1.000000,(j:1.000000,k:1.000000,l:1.000000):1.000000,(d:1.000000,e:1.000000,f:1.000000):1.000000):1.000000);";

  /* and then feed this string to the parser */
  char** taxname_lookup_table = NULL;
  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  Tree* swap_tree = NULL;

  int e_index;
  int min_num_moved=0;

  int max_branches_boot = ref_tree->nb_taxa*2-2;
  int n = ref_tree->nb_taxa;
  int m = ref_tree->nb_edges;
  int i;
  short unsigned** c_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of complements */
  for (i=0; i<m; i++) c_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** i_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of intersections */
  for (i=0; i<m; i++) i_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** hamming = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of Hamming distances */
  for (i=0; i<m; i++) hamming[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned* min_dist = (short unsigned*) malloc(m*sizeof(short unsigned)); /* array of min Hamming distances */

  swap_tree = complete_parse_nh(swap_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  for (i = 0; i < m; i++) {
    min_dist[i] = n; /* initialization to the nb of taxa */
  }
  min_num_moved = 3;
  e_index=8;

  /* calculation of the C and I matrices (see Brehelin/Gascuel/Martin) */
  update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
  update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);
  
  if(min_dist[e_index] > min_num_moved){
    fprintf(stderr,"TRANSFER Test 1 : Error : The min_dist of the swaped branch is > the number of swaped taxa %d>%d\n",min_dist[e_index],min_num_moved);
    exit(EXIT_FAILURE);
  }

  free_tree(swap_tree);
  
  for (i=0; i<m; i++) {
    free(c_matrix[i]);
    free(i_matrix[i]);
    free(hamming[i]);
  }
  free(c_matrix);
  free(i_matrix);
  free(hamming);
  free(min_dist);
  free_tree(ref_tree);

  fprintf(stderr,"TRANSFER Test 1 : OK\n");

  return EXIT_SUCCESS;
}

/**
   We test the TRANSFER Support for branches of a huge multifurcated tree (ncbitax)
   For which we swap branches 100 times
 */
int test_transfer_2(){
  srand(time(NULL));
  char *ref_tree_string = "((((((Dasypus_kappleri,Dasypus_sp.,Dasypus_novemcinctus),Cabassous_unicinctus),(Tamandua_tetradactyla,((Bradypus_tridactylus,Bradypus_torquatus),Choloepus_didactylus))),((Procavia_capensis,Dendrohyrax_dorsalis),Echinops_telfairi,(Elephantulus_sp._VB001,Rhynchocyon_petersi,Macroscelides_proboscideus),(Dugong_dugon,Trichechus_manatus),((Loxodonta_africana,Loxodonta_africana_knochenhaueri,Loxodonta_cyclotis),(Mammuthus_columbi,Mammuthus_primigenius),Elephas_maximus,Elephas_maximus_indicus,Mammut_americanum),Orycteropus_afer,(Eremitalpa_granti,Chrysochloris_asiatica)),((Galeopterus_variegatus,(((Tarsius_syrichta,Tarsius_wallacei,Tarsius_bancanus,Tarsius_dentatus,Tarsius_lariang),((((Nomascus_leucogenys,(Hylobates_agilis,Hylobates_moloch,Hylobates_lar),Symphalangus_syndactylus),(((Homo_sapiens,Homo_sapiens_ssp_Denisova,Homo_heidelbergensis),Gorilla_gorilla,Gorilla_gorilla_gorilla,(Pan_paniscus,Pan_troglodytes,(Pan_troglodytes_troglodytes,Pan_troglodytes_ellioti))),(Pongo_pygmaeus,Pongo_abelii))),(((Colobus_satanas,Colobus_guereza),(Rhinopithecus_brelichi,Rhinopithecus_avunculus,Rhinopithecus_roxellana,Rhinopithecus_bieti_2_RL2012),Semnopithecus_entellus,Nasalis_larvatus,Presbytis_melalophos,Piliocolobus_badius,(Pygathrix_nemaeus,Pygathrix_nigripes),Simias_concolor,(Trachypithecus_francoisi,Trachypithecus_obscurus,Trachypithecus_pileatus,Trachypithecus_johnii,Trachypithecus_cristatus),Procolobus_verus),((Mandrillus_sphinx,Mandrillus_leucophaeus),(Chlorocebus_aethiops,Chlorocebus_sabaeus,Chlorocebus_pygerythrus,Chlorocebus_cynosuros,Chlorocebus_tantalus),Rungwecebus_kipunji,(Lophocebus_aterrimus,Lophocebus_albigena),Allenopithecus_nigroviridis,(Papio_kindae,Papio_hamadryas,Papio_anubis,Papio_ursinus,Papio_papio,Papio_cynocephalus),Erythrocebus_patas,(Miopithecus_talapoin,Miopithecus_ogouensis),(Cercopithecus_dryas,Cercopithecus_aethiops,Cercopithecus_cephus,(Cercopithecus_cephus_cephus,Cercopithecus_cephus_ngottoensis),Cercopithecus_kandti,Cercopithecus_campbelli,Cercopithecus_nictitans,(Cercopithecus_nictitans_nictitans,Cercopithecus_nictitans_martini),(Cercopithecus_ascanius_katangae,Cercopithecus_ascanius_schmidti,Cercopithecus_ascanius_whitesidei),Cercopithecus_roloway,Cercopithecus_erythrotis_camerunensis,Cercopithecus_erythrogaster,Cercopithecus_erythrogaster_pococki,Cercopithecus_petaurista,(Cercopithecus_petaurista_petaurista,Cercopithecus_petaurista_buettikoferi),(Cercopithecus_preussi_insularis,Cercopithecus_preussi_preussi),Cercopithecus_doggetti,Cercopithecus_diana,Cercopithecus_neglectus,Cercopithecus_lhoesti,Cercopithecus_pogonias,(Cercopithecus_pogonias_schwarzianus,Cercopithecus_pogonias_nigripes,Cercopithecus_pogonias_grayi),(Cercopithecus_wolfi_elegans,Cercopithecus_wolfi_pyrogaster),Cercopithecus_albogularis,(Cercopithecus_albogularis_kolbi,Cercopithecus_albogularis_erythrarchus,Cercopithecus_albogularis_monoides,Cercopithecus_albogularis_labiatus,Cercopithecus_albogularis_albotorquatus,Cercopithecus_albogularis_moloneyi,Cercopithecus_albogularis_francescae),Cercopithecus_mitis,(Cercopithecus_mitis_mitis,Cercopithecus_mitis_opisthostictus,Cercopithecus_mitis_boutourlinii,Cercopithecus_mitis_stuhlmanni,Cercopithecus_mitis_heymansi),Cercopithecus_hamlyni,Cercopithecus_solatus,Cercopithecus_mona),(Macaca_nemestrina,Macaca_fascicularis,Macaca_fuscata,Macaca_silenus,Macaca_thibetana,Macaca_assamensis,Macaca_arctoides,Macaca_mulatta,Macaca_tonkeana,Macaca_nigra,Macaca_sylvanus),Theropithecus_gelada,(Cercocebus_atys,Cercocebus_agilis,Cercocebus_torquatus,Cercocebus_chrysogaster)))),((((Chiropotes_israelita,Chiropotes_albinasus),Cacajao_calvus,Pithecia_pithecia),(Callicebus_cupreus,Callicebus_donacophilus,Callicebus_lugens)),(Aotus_nancymaae,Aotus_trivirgatus,Aotus_azarae,Aotus_azarai,Aotus_azarae_azarai,Aotus_lemurinus),((Saimiri_boliviensis,Saimiri_oerstedii,Saimiri_oerstedii_citrinellus,Saimiri_sciureus,Saimiri_sciureus_macrodon),(Callimico_goeldii,Leontopithecus_rosalia,Saguinus_oedipus,(Callithrix_jacchus,Callithrix_pygmaea)),((Cebus_apella,Sapajus_xanthosternos),Cebus_albifrons)),(((Ateles_geoffroyi,Ateles_belzebuth,Ateles_paniscus),Lagothrix_lagotricha,Brachyteles_arachnoides),Alouatta_caraya)))),(Daubentonia_madagascariensis,(((Lepilemur_hubbardorum,Lepilemur_ruficaudatus),Megaladapis_edwardsi),Cheirogaleus_medius,(Indri_indri,Avahi_laniger,(Propithecus_verreauxi,Propithecus_coquereli)),(Lemur_catta,Hapalemur_griseus,Prolemur_simus,(Eulemur_fulvus,Eulemur_rufus,Eulemur_macaco,Eulemur_mongoz,Eulemur_rubriventer),(Varecia_rubra,Varecia_variegata)),Palaeopropithecus_ingens),(((Loris_lydekkerianus,Loris_tardigradus),(Nycticebus_bengalensis,Nycticebus_pygmaeus,Nycticebus_coucang),Perodicticus_potto,Perodicticus_potto_edwarsi),((Otolemur_garnettii,Otolemur_crassicaudatus),(Galago_senegalensis,Galago_moholi),Galagoides_demidoff)))),((((Lepus_comus,Lepus_capensis,Lepus_granatensis,Lepus_arcticus,Lepus_peguensis,Lepus_yarkandensis),Oryctolagus_cuniculus,Sylvilagus_floridanus),(Ochotona_collaris,Ochotona_pallasi,Ochotona_turuchanensis,Ochotona_rufescens,Ochotona_cansus,Ochotona_curzoniae,Ochotona_mantchurica,Ochotona_princeps,Ochotona_hyperborea,Ochotona_sp._WL127_2006,Ochotona_pusilla,Ochotona_dauurica)),(((Anomalurus_pelii,Anomalurus_sp._GP-2005),((Dipodomys_merriami,Dipodomys_ordii,Dipodomys_phillipsii,Dipodomys_microps,Dipodomys_panamintinus,Dipodomys_nelsoni,Dipodomys_deserti,Dipodomys_compactus,Dipodomys_heermanni),((Liomys_irroratus,Liomys_spectabilis,Liomys_salvini,Liomys_pictus),(Heteromys_oresterus,Heteromys_gaumeri,Heteromys_desmarestianus)),((Chaetodipus_penicillatus,Chaetodipus_hispidus,Chaetodipus_baileyi,Chaetodipus_californicus,Chaetodipus_eremicus,Chaetodipus_arenarius,Chaetodipus_pernix,Chaetodipus_intermedius,Chaetodipus_spinatus,Chaetodipus_formosus),(Perognathus_parvus,Perognathus_flavescens,Perognathus_flavus,Perognathus_longimembris,Perognathus_merriami))),(Muscardinus_avellanarius,Glis_glis),((Orthogeomys_heterodus,Orthogeomys_grandis,Orthogeomys_hispidus,Orthogeomys_matagalpae,Orthogeomys_underwoodi,Orthogeomys_cavator),Pappogeomys_bulleri,Geomys_breviceps,Zygogeomys_trichopus,(Thomomys_umbrinus,Thomomys_talpoides,Thomomys_sheldoni,Thomomys_bulbivorus),(Cratogeomys_perotensis,Cratogeomys_neglectus,Cratogeomys_castanops,Cratogeomys_gymnurus,Cratogeomys_fulvescens,Cratogeomys_fumosus,Cratogeomys_tylorhinus,Cratogeomys_zinseri,Cratogeomys_goldmani)),(Castor_fiber,Castor_canadensis),(((Zapus_princeps,Zapus_trinotatus,Zapus_hudsonius),Eozapus_setchuanus,Napaeozapus_insignis),(Jaculus_jaculus,Dipus_sagitta),Sicista_concolor,(Allactaga_sibirica,Allactaga_toussi,Allactaga_elater,Allactaga_firouzi),Euchoreutes_naso),(((Spermophilopsis_leptodactylus,Xerus_erythropus),(Spermophilus_parryii,Spermophilus_lateralis,(Cynomys_leucurus,Cynomys_ludovicianus),(Spermophilus_erythrogenys,Spermophilus_alashanicus,Spermophilus_major,Spermophilus_pygmaeus,Spermophilus_suslicus,Spermophilus_dauricus),Spermophilus_tridecemlineatus,(Tamias_minimus,Tamias_sibiricus,Tamias_striatus,Tamias_amoenus),(Marmota_himalayana,Marmota_monax)),Heliosciurus_gambianus),(((Hylopetes_spadiceus,Hylopetes_phayrei,Hylopetes_alboniger),Pteromys_volans,Petinomys_setosus,Glaucomys_volans),Tamiasciurus_hudsonicus,(Sciurus_carolinensis,Sciurus_vulgaris)),((Funambulus_palmarum,Funisciurus_anerythrus),((Callosciurus_sp._1_MG2013,Callosciurus_notatus,Callosciurus_erythraeus),Exilisciurus_exilis,Tamiops_swinhoei,Dremomys_rufigenis)),Ratufa_bicolor),((Acomys_cahirinus,((Maxomys_surifer,Maxomys_moi,Maxomys_rajah,Maxomys_whiteheadi),Millardia_meltada,Berylmys_bowersi,(Mastomys_erythroleucus,Mastomys_natalensis,Mastomys_kollmannspergeri,Mastomys_coucha,Mastomys_huberti),((Mus_musculus,Mus_musculus_domesticus,Mus_spretus,Mus_fragilicauda),Mus_pahari),Leggadina_lakedownensis,(Malacomys_cansdalei,Malacomys_longipes),Chiromyscus_chiropus,Heimyscus_fumosus,(Myomyscus_verreauxii,Myomyscus_brockmani),(Hylomyscus_kaimosae,Hylomyscus_grandis,Hylomyscus_walterverheyeni,Hylomyscus_pamfi,(Hylomyscus_sp._1,Hylomyscus_sp._6,Hylomyscus_sp._2),Hylomyscus_aeta,Hylomyscus_alleni,Hylomyscus_stella,Hylomyscus_simus,Hylomyscus_parvus),Pseudomys_chapmani,(Rhabdomys_pumilio,Rhabdomys_dilectus),(Praomys_hartwigi,Praomys_misonnei,Praomys_sp._A,Praomys_morio,Praomys_delectorum,Praomys_tullbergi,Praomys_derooi,Praomys_rostratus,Praomys_daltoni),(Apodemus_peninsulae,Apodemus_draco,Apodemus_flavicollis,Apodemus_latronum,Apodemus_agrarius,Apodemus_chejuensis),Bandicota_indica,(Leopoldamys_neilli,Leopoldamys_sabanus,Leopoldamys_edwardsi),(Micromys_erythrotis,Micromys_minutus),(Rattus_leucopus,Rattus_exulans,Rattus_norvegicus,(Rattus_sp._ABTC_42808,Rattus_sp._ABTC_47998,Rattus_sp._abtc_43216,Rattus_sp._abtc_45409),Rattus_rattus,Rattus_losea,Rattus_tanezumi,Rattus_tanezumi_sladeni,Rattus_giluwensis,Rattus_niobe,Rattus_fuscipes,Rattus_argentiventer,Rattus_andamanensis,Rattus_remotus,Rattus_tiomanicus),Niviventer_confucianus),(Tatera_indica,Rhombomys_opimus,(Gerbillus_nanus,Gerbillus_sp._1_TCB2013),(Meriones_unguiculatus,Meriones_libycus,Meriones_meridianus))),Typhlomys_cinereus,(((Scotinomys_xerampelinus,Scotinomys_teguina),Habromys_lophurus,Neotomodon_alstoni,(Peromyscus_melanocarpus,Peromyscus_mayensis,Peromyscus_levipes,Peromyscus_leucopus,Peromyscus_stirtoni,Peromyscus_grandis,Peromyscus_maniculatus,Peromyscus_pectoralis,Peromyscus_aztecus,Peromyscus_mexicanus,Peromyscus_yucatanicus),Osgoodomys_banderanus,(Reithrodontomys_fulvescens,Reithrodontomys_microdon,Reithrodontomys_spectabilis,Reithrodontomys_sumichrasti,Reithrodontomys_mexicanus,Reithrodontomys_gracilis,Reithrodontomys_megalotis),Neotoma_cinerea,Isthmomys_pirrensis),(Lasiopodomys_mandarinus,(Eothenomys_melanogaster,Eothenomys_chinensis,Eothenomys_proditor,Eothenomys_custos,Eothenomys_sp._2_SL-2010a,Eothenomys_inez,Eothenomys_eleusis,Eothenomys_miletus,Eothenomys_eva),(Neodon_irene,Neodon_leucurus),Lemmus_trimucronatus,Phenacomys_intermedius,(Myodes_glareolus,Myodes_gapperi,Myodes_rufocanus,Myodes_rutilus),(Microtus_pennsylvanicus,Microtus_middendorffii,Microtus_guatemalensis,Microtus_ochrogaster,Microtus_longicaudus,Microtus_limnophilus,Microtus_kikuchii,Microtus_rossiaemeridionalis,Microtus_levis),(Dicrostonyx_groenlandicus,Dicrostonyx_richardsoni),Arvicola_terrestris,Alticola_stracheyi),(Thaptomys_nigrita,Rheomys_thomasi,Wiedomys_cerradensis,Sooretamys_angouya,(Oecomys_bicolor,Oecomys_cf._rex,Oecomys_sp._CMV2014,Oecomys_auyantepui),(Bolomys_lasiurus,Necromys_urichi),(Delomys_dorsalis,Delomys_sublineatus),(Oryzomys_melanotis,Oryzomys_alfaroi),(Oryzomys_yunganus,Oryzomys_capito,Hylaeamys_megacephalus,Oryzomys_megacephalus,Oryzomys_perenensis),Deltamys_kempi,Brucepattersonius_soricinus,Sigmodon_hispidus,Scolomys_melanops,Zygodontomys_brevicauda,Oryzomys_couesi,Nectomys_squamipes,(Neacomys_sp.,Neacomys_spinosus,Neacomys_guianae,Neacomys_paracou),(Rhipidomys_macconnelli,Rhipidomys_leucodactylus,Rhipidomys_nitela),Juliomys_pictipes,(Euryoryzomys_macconnelli,Oryzomys_macconnelli,Euryoryzomys_russatus),(Oligoryzomys_fulvescens,Oligoryzomys_nigripes,Oligoryzomys_fornesi,Oligoryzomys_flavescens),(Akodon_serrensis,Akodon_cursor,Akodon_montensis,Akodon_azarae),Euneomys_mordax,Calomys_expulsus,Oryzomys_albigularis),(Tylomys_nudicaudus,Ototylomys_phyllotis),(Tscherskia_triton,Mesocricetus_auratus,Allocricetulus_curtatus,(Cricetulus_griseus,Cricetulus_migratorius,Cricetulus_longicaudatus,Cricetulus_kamensis),(Phodopus_campbelli,Phodopus_roborovskii))),(Saccostomus_campestris,(Cricetomys_gambianus,Cricetomys_emini)),((Rhizomys_pruinosus,Rhizomys_sinensis),(Nannospalax_golani,Nannospalax_galili,Nannospalax_ehrenbergi,Nannospalax_judaei),((Myospalax_aspalax,Myospalax_psilurus),(Eospalax_baileyi,Eospalax_cansus,Eospalax_rothschildi))))),(Hydrochoerus_hydrochaeris,(Ctenomys_pearsoni,Ctenomys_lami,Ctenomys_dorbignyi,Ctenomys_perrensi,Ctenomys_torquatus,Ctenomys_sociabilis,Ctenomys_minutus,Ctenomys_conoveri,Ctenomys_rionegrensis,Ctenomys_leucodon),Cavia_porcellus,Heterocephalus_glaber,Cuniculus_paca,Chinchilla_lanigera,Thryonomys_swinderianus,Hystrix_indica,(Makalata_didelphoides,(Phyllomys_dasythrix,Phyllomys_blainvillii,Phyllomys_pattoni,Phyllomys_nigrispinus,Phyllomys_sp._ACL-2011),(Proechimys_hoplomyoides,Proechimys_cuvieri,Proechimys_simonsi,Proechimys_sp._bkl1,Proechimys_quadruplicatus,Proechimys_longicaudatus,Proechimys_guyannensis,Proechimys_gularis),Euryzygomatomys_spinosus,Echimys_semivillosus,Mesomys_hispidus,Trinomys_dimidiatus),Dasyprocta_leporina,(Spalacopus_cyanus,Octomys_mimax,Tympanoctomys_barrerae,Octodon_degus),(Erethizon_dorsata,Coendou_insidiosus)))),(Tupaia_minor,Tupaia_belangeri)),(((Tapirus_indicus,Tapirus_terrestris),((Rhinoceros_unicornis,Rhinoceros_sondaicus),Dicerorhinus_sumatrensis,Diceros_bicornis,Coelodonta_antiquitatis,Ceratotherium_simum),((Equus_burchellii,Equus_zebra),Equus_grevyi,(Equus_asinus_somalicus,Equus_asinus_africanus),(Equus_ferus_caballus,Equus_caballus,Equus_ferus_przewalskii))),((Uropsilus_soricipes,Talpa_europaea,Scapanulus_oweni,Condylura_cristata,Mogera_wogura,Galemys_pyrenaicus,Neurotrichus_gibbsii,Urotrichus_talpoides),((Neotetracus_sinensis,Echinosorex_gymnura,Hylomys_suillus,Neohylomys_hainanensis),(Erinaceus_europaeus,Hemiechinus_auritus)),(((Suncus_megalura,Suncus_murinus),(Crocidura_flavescens,Crocidura_muricauda,Crocidura_olivieri,Crocidura_fuliginosa,Crocidura_brunnea,Crocidura_buettikoferi,Crocidura_attenuata,Crocidura_shantungensis,Crocidura_grandiceps,Crocidura_cf._tanakae,Crocidura_douceti,Crocidura_viaria,Crocidura_jouvenetae,Crocidura_wuchihensis,Crocidura_obscurior,Crocidura_nimbasilvanus,Crocidura_goliath_nimbasilvanus,Crocidura_russula)),((Episoriculus_caudatus,Episoriculus_fumidus),Neomys_fodiens,Nectogale_elegans,Anourosorex_squamipes,(Blarinella_griselda,Blarinella_quadraticauda),(Chodsigoa_parca,Soriculus_sodalis),(Sorex_tundrensis,Sorex_isodon,Sorex_minutissimus,Sorex_cylindricauda,Sorex_bedfordiae,Sorex_arcticus,Sorex_cinereus,Sorex_caecutiens,Sorex_unguiculatus,Sorex_trowbridgii,Sorex_fumeus),Blarina_brevicauda,Soriculus_nigrescens))),(Manis_pentadactyla,Smutsia_gigantea,Manis_javanica,Phataginus_tricuspis,Phataginus_tetradactyla,Smutsia_temminckii),(((Giraffa_camelopardalis,Giraffa_camelopardalis_rothschildi,(Moschus_moschiferus,Moschus_anhuiensis,Moschus_berezovskii),(Hydropotes_inermis,((Mazama_sp.,Mazama_nemorivaga),(Alces_americanus,Alces_alces),(Odocoileus_hemionus,Odocoileus_virginianus),Capreolus_capreolus,Rangifer_tarandus),((Muntiacus_crinifrons,Muntiacus_reevesi,Muntiacus_muntjak),Elaphodus_cephalophus),((Cervus_unicolor,Rusa_unicolor,Rusa_alfredi,Rusa_timorensis),Elaphurus_davidianus,Przewalskium_albirostris,(Dama_mesopotamica,Dama_dama),(Cervus_canadensis,Cervus_nippon,(Cervus_nippon_yakushimae,Cervus_hortulorum,Cervus_yesoensis),Cervus_elaphus),Rucervus_duvaucelii,(Axis_axis,Axis_porcinus))),(Redunca_redunca,((Hemitragus_jemlahicus,Hemitragus_jayakari),Rupicapra_rupicapra,Oreamnos_americanus,Budorcas_taxicolor,(Capra_sibirica,Capra_falconeri,Capra_hircus,Capra_nubiana),Ammotragus_lervia,(Ovis_dalli,Ovis_aries,Ovis_canadensis,Ovis_ammon_hodgsoni),Ovibos_moschatus,(Capricornis_crispus,Capricornis_milneedwardsii),(Naemorhedus_griseus,Naemorhedus_caudatus,Naemorhedus_goral,Naemorhedus_baileyi)),Aepyceros_melampus,(Antidorcas_marsupialis,Nanger_granti,Saiga_tatarica,Pantholops_hodgsonii,Antilope_cervicapra,Ourebia_ourebi,(Neotragus_pygmaeus,Neotragus_batesi,Neotragus_moschatus),(Gazella_erlangeri,Gazella_subgutturosa),Procapra_gutturosa,Eudorcas_rufifrons),(Beatragus_hunteri,(Damaliscus_lunatus,Damaliscus_pygargus,Damaliscus_pygargus_phillipsi)),(Hippotragus_equinus,(Oryx_dammah,Oryx_gazella)),(Tetracerus_quadricornis,(Bubalus_bubalis,Bubalus_arnee,Bubalus_carabanensis),(Bos_taurus_indicus,Bos_gaurus,Bos_taurus),Pseudoryx_nghetinhensis,Bison_bonasus,(Tragelaphus_angasii,Tragelaphus_oryx,Tragelaphus_imberbis,Tragelaphus_eurycerus,Tragelaphus_eurycerus_eurycerus,Tragelaphus_scriptus),Syncerus_caffer),((Cephalophus_dorsalis,Cephalophus_weynsi,Cephalophus_niger,Cephalophus_ogilbyi,Cephalophus_silvicultor,Cephalophus_natalensis,Cephalophus_nigrifrons,Cephalophus_harveyi,Cephalophus_jentinki,Cephalophus_zebra,Cephalophus_callipygus,Cephalophus_adersi),(Cephalophus_monticola,Philantomba_monticola,Philantomba_maxwellii))),Antilocapra_americana),(Hyemoschus_aquaticus,Moschiola_indica)),((Vicugna_vicugna,Vicugna_pacos,Lama_pacos),(Lama_glama,Lama_guanicoe),(Camelus_bactrianus,Camelus_dromedarius)),(Hippopotamus_amphibius,Hexaprotodon_liberiensis),(((Megaptera_novaeangliae,(Balaenoptera_borealis,Balaenoptera_omurai,Balaenoptera_acutorostrata,Balaenoptera_physalus,Balaenoptera_bonaerensis,Balaenoptera_musculus,Balaenoptera_edeni)),Eschrichtius_robustus,Caperea_marginata,Balaena_mysticetus),((Monodon_monoceros,Delphinapterus_leucas),(Physeter_catodon,Physeter_macrocephalus,Kogia_breviceps),Pontoporia_blainvillei,Platanista_minor,(Phocoenoides_dalli,(Phocoena_phocoena,Phocoena_sinus,Phocoena_spinipinnis),(Neophocaena_asiaeorientalis,Neophocaena_phocaenoides)),(Inia_araguaiaensis,Inia_boliviensis,Inia_geoffrensis),((Mesoplodon_stejnegeri,Mesoplodon_densirostris),Hyperoodon_ampullatus,Ziphius_cavirostris,Berardius_bairdii),(Orcinus_orca,Pseudorca_crassidens,Lissodelphis_borealis,Grampus_griseus,(Stenella_frontalis,Stenella_attenuata,Stenella_coeruleoalba),(Orcaella_heinsohni,Orcaella_brevirostris),(Tursiops_australis,Tursiops_aduncus,Tursiops_truncatus),Cephalorhynchus_heavisidii,(Globicephala_melas,Globicephala_macrorhynchus),Peponocephala_electra,Sotalia_fluviatilis,Steno_bredanensis,(Lagenorhynchus_obliquidens,Lagenorhynchus_albirostris,Lagenorhynchus_acutus),Sousa_chinensis),Lipotes_vexillifer)),((Tayassu_pecari,Tayassu_tajacu,Pecari_tajacu),(Potamochoerus_porcus,(Sus_verrucosus,Sus_cebifrons,Sus_scrofa,Sus_scrofa_domesticus)))),((Sphaerias_blanfordi,(Eonycteris_spelaea,(Macroglossus_sobrinus,Macroglossus_minimus),(Melonycteris_woodfordi,Melonycteris_fardoulisi,Melonycteris_melanops)),((Megaloglossus_woermanni,(Micropteropus_pusillus,Epomophorus_labiatus)),Dobsonia_inermis,(Megaerops_ecaudatus,Megaerops_niphanae),Eidolon_helvum,Balionycteris_maculata,Dyacopterus_spadiceus,(Pteropus_samoensis,Pteropus_dasymallus,Pteropus_giganteus,Pteropus_vampyrus),(Cynopterus_JLE_sp._A,Cynopterus_brachyotis,Cynopterus_sphinx,Cynopterus_titthaecheilus,Cynopterus_horsfieldii),(Rousettus_aegyptiacus,Rousettus_amplexicaudatus,Rousettus_leschenaultii),Chironax_melanocephalus)),((Brachyphylla_cavernarum,((Glyphonycteris_daviesi,Glyphonycteris_sylvestris),Macrophyllum_macrophyllum,(Lophostoma_carrikeri,Lophostoma_brasiliense,Lophostoma_schulzi,Lophostoma_silvicolum),(Trachops_cirrhosus_PS3,Trachops_cirrhosus_PS1,Trachops_cirrhosus),(Lonchorhina_aurita,Lonchorhina_inusitata),(Phylloderma_stenops_PS1,Phylloderma_stenops),(Phyllostomus_latifolius,Phyllostomus_elongatus,Phyllostomus_hastatus,Phyllostomus_discolor),Tonatia_saurophila,Vampyrum_spectrum,(Trinycteris_nicefori,Micronycteris_hirsuta,Micronycteris_schmidtorum,Lampronycteris_brachyotis,Micronycteris_microtis,Micronycteris_megalotis,Micronycteris_minuta),Chrotopterus_auritus,(Mimon_cozumelae,Mimon_crenulatum)),(Mesophylla_macconnelli,(Platyrrhinus_helleri,Platyrrhinus_infuscus,Platyrrhinus_vittatus),(Sturnira_magna,Sturnira_erythromos,Sturnira_tildae,Sturnira_lilium,Sturnira_bidens,Sturnira_luisi,Sturnira_ludovici),(Vampyressa_brocki,Vampyressa_bidens,Vampyressa_thyone),(Artibeus_intermedius,Artibeus_aztecus,Artibeus_obscurus,Artibeus_cinereus,Artibeus_phaeotis,Artibeus_jamaicensis,(Artibeus_watsoni,Artibeus_gnomus),Artibeus_planirostris,Artibeus_lituratus,Artibeus_concolor,Artibeus_bogotensis,Artibeus_anderseni),(Chiroderma_villosum,Chiroderma_doriae,Chiroderma_trinitatum),Uroderma_bilobatum,Ametrida_centurio,Vampyrodes_caraccioli),((Carollia_brevicauda,Carollia_sowelli,Carollia_perspicillata,Carollia_brevicauda_PS1,Carollia_castanea),Rhinophylla_pumilio),(Lonchophylla_thomasi,Lionycteris_spurrelli,(Lonchophylla_robusta,Lonchophylla_chocoana)),(Hylonycteris_underwoodi,(Choeroniscus_sp.,Choeroniscus_minor,Choeroniscus_godmani),(Glossophaga_longirostris,Glossophaga_commissarisi,Glossophaga_soricina),(Anoura_geoffroyi,Anoura_caudifer,Anoura_latidens)),(Diaemus_youngi,Diphylla_ecaudata,Desmodus_rotundus)),Furipterus_horrens,(Nycteris_tragata,Nycteris_thebaica),Mystacina_tuberculata,((Pteronotus_parnellii,Pteronotus_personatus,Pteronotus_rubiginosus,Pteronotus_gymnonotus),Mormoops_megalophylla),(Rhogeessa_io,Hesperoptenus_tickelli,(Nyctalus_lasiopterus,Nyctalus_leisleri),Philetor_brachypterus,(Barbastella_barbastellus,Barbastella_leucomelas,Barbastella_darjelingensis),Parastrellus_hesperus,(Scotophilus_dinganii,Scotophilus_kuhlii,Scotophilus_heathii),Glischropus_tylopus,Nycticeius_humeralis,Chalinolobus_tuberculatus,(Miniopterus_magnater,Miniopterus_fuliginosus,Miniopterus_pusillus,Miniopterus_medius),Harpiola_isodon,Ia_io,(Antrozous_pallidus,Bauerus_dubiaquercus),Corynorhinus_townsendii,Euderma_maculatum,(Kerivoula_cf._papillosa,Kerivoula_titania,Kerivoula_hardwickii,Kerivoula_pellucida,Kerivoula_cf._lenis,Kerivoula_minuta,Kerivoula_sp._FAK-2010,Kerivoula_kachinensis,Kerivoula_picta,Kerivoula_papillosa,Kerivoula_intermedia,Kerivoula_cf._hardwickii),(Myotis_yumanensis,Myotis_ciliolabrum,Myotis_formosus,Myotis_riparius,Myotis_hasseltii,Myotis_macrotarsus,Myotis_volans,Myotis_aurascens,Myotis_cf._aurascens,Myotis_californicus,Myotis_blythii,Myotis_blythii_omari,Myotis_petax,Myotis_frater,Myotis_septentrionalis,Myotis_cf._alcathoe,Myotis_daubentonii,Myotis_phanluongi,Myotis_brandtii,Myotis_bombinus,Myotis_bechsteinii,Myotis_capaccinii,Myotis_cf._laniger,Myotis_riparius_PS3,Myotis_schaubi,Myotis_rosseti,Myotis_nigricans_PS2,Myotis_pilosus,Myotis_gomantongensis,Myotis_montivagus,Myotis_siligorensis,Myotis_annamiticus,Myotis_cf._muricola,Myotis_keaysi,Myotis_davidii,Myotis_myotis,Myotis_evotis,Myotis_ikonnikovi,Myotis_riparius_PS2,Myotis_annectans,Myotis_muricola,Myotis_dasycneme,Myotis_horsfieldii,Myotis_macrodactylus,Myotis_nattereri,Myotis_laniger,Myotis_lucifugus,Myotis_taiwanensis,Myotis_mystacinus,Myotis_sodalis,Myotis_chinensis,Myotis_annatessae,Myotis_alcathoe,Myotis_nigricans_PS1,Myotis_albescens),Lasionycteris_noctivagans,(Hypsugo_cadornae,Hypsugo_crassulus_bellieri,Pipistrellus_eisentrauti,Hypsugo_pulveratus),(Lasiurus_atratus,Lasiurus_intermedius,Lasiurus_blossevillii,Lasiurus_borealis,Lasiurus_cinereus,Lasiurus_seminolus,Lasiurus_xanthinus),(Murina_aenea,Murina_fionae,Murina_ussuriensis,Murina_harpioloides,Murina_leucogaster,Murina_hilgendorfi,Murina_tubinaris,Murina_sp.,Murina_lorelieae,Murina_cyclotis,Murina_walstoni,Murina_cf._cyclotis,Murina_huttoni,Murina_eleryi,Murina_annamitica),(Eptesicus_furinalis,Eptesicus_JLE_sp._A,Eptesicus_nasutus,Eptesicus_nilssonii,Eptesicus_bottae_anatolicus,Eptesicus_fuscus,Eptesicus_chiriquinus,Eptesicus_serotinus),(Plecotus_macrobullaris,Plecotus_macrobullaris_alpinus,Plecotus_auritus,Plecotus_austriacus,Plecotus_cf._strelkovi,Plecotus_rafinesquii,Plecotus_kolombatovici,Plecotus_teneriffae_gaisleri,Plecotus_christii,Plecotus_ognevi),(Tylonycteris_robustula,Tylonycteris_pachypus),Otonycteris_hemprichii,Harpiocephalus_harpia,(Pipistrellus_pygmaeus,Pipistrellus_javanicus,Pipistrellus_cf._coromandra,Pipistrellus_abramus,Pipistrellus_kuhlii,Pipistrellus_kuhlii_kuhlii,Pipistrellus_sp._MIBZPL02288,Pipistrellus_nathusii,Pipistrellus_coromandra,Pipistrellus_kuhlii_deserti,Pipistrellus_deserti,Pipistrellus_paterculus,Pipistrellus_subflavus,Pipistrellus_rueppellii,Pipistrellus_pipistrellus,Pipistrellus_tenuis),Eudiscopus_denticulus,(Neoromicia_capensis,Eptesicus_brunneus),Vespertilio_murinus,Scotomanes_ornatus,Idionycteris_phyllotis),(Chaerephon_nigeriae,Cynomops_paranus,Nyctinomops_laticaudatus,Tadarida_brasiliensis,Eumops_hansae,(Molossus_molossus,Molossus_rufus),Molossops_neglectus),(Megaderma_spasma,Megaderma_lyra),(Asellia_tridens,(Hipposideros_ruber,Hipposideros_cineraceus,(Hipposideros_larvatus,Hipposideros_grandis,Hipposideros_cf._larvatus),Hipposideros_khaokhouayensis,Hipposideros_cf._bicolor,Hipposideros_bicolor_131,Hipposideros_bicolor31,Hipposideros_ater,Hipposideros_ridleyi,Hipposideros_cyclops,Hipposideros_beatus,Hipposideros_cervinus,Hipposideros_diadema,Hipposideros_galeritus,Hipposideros_lylei,Hipposideros_hypophyllus,Hipposideros_commersoni,Hipposideros_speoris,Hipposideros_armiger,Hipposideros_pratti,Hipposideros_pomona,Hipposideros_rotalis,Hipposideros_CMF_sp._C),Aselliscus_stoliczkanus,Coelops_frithii),(Noctilio_albiventris,Noctilio_leporinus,Noctilio_albiventris_PS2),(Rhinolophus_sinicus,Rhinolophus_rex,Rhinolophus_marshalli,Rhinolophus_stheno,Rhinolophus_creaghi,Rhinolophus_lepidus,Rhinolophus_luctus,Rhinolophus_alcyone,Rhinolophus_shameli,Rhinolophus_borneensis,Rhinolophus_affinis,Rhinolophus_malayanus,Rhinolophus_beddomei,Rhinolophus_cognatus,Rhinolophus_paradoxolophus,Rhinolophus_formosae,Rhinolophus_rouxii,Rhinolophus_yunnanensis,Rhinolophus_hipposideros,Rhinolophus_pearsonii,Rhinolophus_yunanensis,Rhinolophus_hildebrandtii,Rhinolophus_cf._lepidus,Rhinolophus_chaseni,Rhinolophus_pusillus,Rhinolophus_philippinensis,Rhinolophus_trifoliatus,Rhinolophus_clivosus,Rhinolophus_coelophyllus,Rhinolophus_euryale,Rhinolophus_cf._pusillus,Rhinolophus_cf._thomasi,Rhinolophus_macrotis,Rhinolophus_acuminatus,Rhinolophus_ferrumequinum,(Rhinolophus_ferrumequinum_quelpartis,Rhinolophus_ferrumequinum_korai)),(((Diclidurus_isabellus,Diclidurus_albus),Cyttarops_alecto,(Emballonura_raffrayana,Emballonura_serii,Emballonura_monticola,Emballonura_beccarii,Emballonura_semicaudata,Emballonura_alecto),Cormura_brevirostris,Mosia_nigrescens,Rhynchonycteris_naso,(Saccopteryx_bilineata,Saccopteryx_leptura,Saccopteryx_canescens),(Balantiopteryx_io,Balantiopteryx_plicata),(Peropteryx_macrotis,Peropteryx_leucoptera)),(Taphozous_longimanus,Taphozous_sp._CS-2014,Taphozous_melanopogon)),(Thyroptera_tricolor,Thyroptera_lavali))),((Nandinia_binotata,((Crossarchus_obscurus,Crossarchus_platycephalus),Ichneumia_albicauda,Liberiictis_kuhni,Herpestes_javanicus,Atilax_paludinosus),((Viverricula_indica,Civettictis_civetta,Genetta_servalina),Paradoxurus_hermaphroditus,Prionodon_pardicolor),(((Lynx_lynx,Lynx_canadensis),Leopardus_wiedii,(Felis_silvestris,Felis_silvestris_silvestris,Felis_catus),(Prionailurus_viverrinus,Prionailurus_planiceps),Puma_concolor),Acinonyx_jubatus,(Panthera_tigris,(Panthera_tigris_amoyensis,Panthera_tigris_corbetti,Panthera_tigris_tigris),Panthera_pardus,Panthera_onca,Panthera_leo,Panthera_leo_persica)),(Crocuta_crocuta,Hyaena_hyaena)),(((Arctocephalus_australis,Arctocephalus_forsteri),Otaria_flavescens,Callorhinus_ursinus,Phocarctos_hookeri,Eumetopias_jubatus),(Spilogale_putorius,Mephitis_mephitis,(Conepatus_chinga,Conepatus_semistriatus)),(Bassariscus_astutus,Potos_flavus,(Nasua_nasua,Nasua_narica),(Procyon_cancrivorus,Procyon_lotor)),(Ailuropoda_melanoleuca,Tremarctos_ornatus,Melursus_ursinus,Arctodus_simus,Helarctos_malayanus,(Ursus_deningeri,Ursus_spelaeus,Ursus_thibetanus,Ursus_thibetanus_mupinensis,Ursus_americanus,Ursus_arctos,Ursus_maritimus)),Odobenus_rosmarus,Odobenus_rosmarus_rosmarus,Ailurus_fulgens,(Cerdocyon_thous,(Lycalopex_culpaeus,Lycalopex_fulvipes,Lycalopex_gymnocercus,Lycalopex_griseus,Lycalopex_vetulus),(Canis_lupus,(Canis_lupus_chanco,Canis_familiaris,Canis_lupus_familiaris),Canis_mesomelas_elongae,Canis_aureus,Canis_latrans,Canis_adustus),(Vulpes_lagopus,Vulpes_corsac,Vulpes_zerda,Vulpes_vulpes,Vulpes_macrotis),Chrysocyon_brachyurus,Urocyon_cinereoargenteus),(Mirounga_angustirostris,Pusa_hispida,Halichoerus_grypus,Ommatophoca_rossii,Phoca_vitulina,Leptonychotes_weddellii,Cystophora_cristata,Erignathus_barbatus),((Arctonyx_collaris,(Meles_meles,Meles_anakuma)),(Neovison_vison,(Mustela_frenata,Mustela_nivalis,Mustela_sibirica,Mustela_putorius,Mustela_erminea,Mustela_kathiah,Mustela_nigripes,Mustela_altaica)),(Martes_pennanti,Martes_flavigula,Martes_americana,Martes_martes),(Hydrictis_maculicollis,Lutra_lutra,(Lontra_canadensis,Lontra_longicaudis),Enhydra_lutris,Pteronura_brasiliensis),Taxidea_taxus,(Galictis_vittata,Galictis_cuja),Melogale_moschata)))))),(Dromiciops_gliroides,Lestoros_inca,(Macrotis_lagotis,Isoodon_macrourus,Perameles_gunnii),((Potorous_gilbertii,Potorous_tridactylus,Potorous_tridactylus_apicalis),Distoechurus_pennatus,Trichosurus_vulpecula,Tarsipes_rostratus,(Burramys_parvus,Cercartetus_nanus),Phascolarctos_cinereus,Pseudocheirus_peregrinus,(Lagorchestes_hirsutus,(Petrogale_lateralis,Petrogale_burbidgei,Petrogale_xanthopus_celeris,Petrogale_rothschildi),Lagostrophus_fasciatus,Macropus_robustus,Thylogale_thetis,(Dendrolagus_lumholtzi,Dendrolagus_goodfellowi)),Vombatus_ursinus,(Dactylopsila_trivirgata,Gymnobelideus_leadbeateri,Petaurus_breviceps)),((Phascogale_tapoatafa,Sarcophilus_harrisii,Dasyuroides_byrnei,(Sminthopsis_crassicaudata,Sminthopsis_douglasi),Parantechinus_apicalis,(Dasyurus_geoffroii,Dasyurus_hallucatus),Antechinus_flavipes),Myrmecobius_fasciatus,Thylacinus_cynocephalus),((Thylamys_elegans,(Micoureus_regina,Micoureus_demerarae,Micoureus_paraguayanus),(Marmosa_murina,Marmosa_waterhousei,Marmosa_mexicana),Metachirus_nudicaudatus,(Philander_opossum,Philander_andersoni),(Didelphis_albiventris,Didelphis_virginiana,Didelphis_aurita,Didelphis_marsupialis),(Marmosops_pinheiroi,Marmosops_parvidens,Marmosops_incanus,Marmosops_noctivagus),Gracilinanus_microtarsus,(Monodelphis_americana,Monodelphis_domestica,Monodelphis_brevicaudata)),Caluromys_philander),Notoryctes_typhlops)),((Tachyglossus_aculeatus,Zaglossus_bruijni),Ornithorhynchus_anatinus));";
//((a:1,b:1,c:1):1,(d:1,e:1,f:1):1,((g:1,h:1,i:1):1,(j:1,k:1,l:1):1,(m:1,n:1,o:1):1):1);"; 

  int r;
  /* and then feed this string to the parser */
  char** taxname_lookup_table = NULL;
  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  Tree* swap_tree = NULL;

  Edge *e;
  int e_index;
  int min_num_moved=0;

  int max_branches_boot = ref_tree->nb_taxa*2-2;
  int n = ref_tree->nb_taxa;
  int m = ref_tree->nb_edges;
  int i;
  short unsigned** c_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of complements */
  for (i=0; i<m; i++) c_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** i_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of intersections */
  for (i=0; i<m; i++) i_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** hamming = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of Hamming distances */
  for (i=0; i<m; i++) hamming[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned* min_dist = (short unsigned*) malloc(m*sizeof(short unsigned)); /* array of min Hamming distances */

  for(r=0;r<100;r++){
    swap_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
    for (i = 0; i < m; i++) {
      min_dist[i] = n; /* initialization to the nb of taxa */
    }

    e = NULL;
    e_index = 0;
    while(e == NULL ||
    	  e->right->nneigh == 1 ||
    	  e->left->nneigh == 1 ){
      e_index = rand_to(swap_tree->nb_edges);
      e = swap_tree->a_edges[e_index];
    }
    min_num_moved = swap_branches(swap_tree,e);

    /* calculation of the C and I matrices (see Brehelin/Gascuel/Martin) */
    update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
    update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);

    if(min_dist[e_index] > min_num_moved){
      fprintf(stderr,"TRANSFER Test 2 after branch swap : Error : The min_dist of the swaped branch is > the number of swaped taxa\n");
      exit(EXIT_FAILURE);
    }
    free_tree(swap_tree);
    swap_tree = NULL;
  }
  
  for (i=0; i<m; i++) {
    free(c_matrix[i]);
    free(i_matrix[i]);
    free(hamming[i]);
  }
  free(c_matrix);
  free(i_matrix);
  free(hamming);
  free(min_dist);

  fprintf(stderr,"TRANSFER Test 2 : OK\n");

  return EXIT_SUCCESS;
}



/**
   We test the TRANSFER Support for branches of the initial tree compared to another tree
 */
int test_transfer_3(){
  srand(time(NULL));
  char *ref_tree_string = "(a:1,b:1,c:1,d:1,e:1,f:1,(g:1,h:1,i:1,j:1,k:1,l:1):1);";
  char *swap_tree_string = "(a:1,b:1,c:1,d:1,h:1,g:1,(f:1,e:1,i:1,j:1,k:1,l:1):1);";
  char *swap_tree2_string = "(a:1,b:1,c:1,i:1,h:1,g:1,(f:1,e:1,d:1,j:1,k:1,l:1):1);";

  /* and then feed this string to the parser */
  char** taxname_lookup_table = NULL;
  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  Tree* swap_tree = NULL;

  int e_index;
  int min_num_moved=0;

  int max_branches_boot = ref_tree->nb_taxa*2-2;
  int n = ref_tree->nb_taxa;
  int m = ref_tree->nb_edges;
  int i;
  short unsigned** c_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of complements */
  for (i=0; i<m; i++) c_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** i_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of intersections */
  for (i=0; i<m; i++) i_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned** hamming = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of Hamming distances */
  for (i=0; i<m; i++) hamming[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
  short unsigned* min_dist = (short unsigned*) malloc(m*sizeof(short unsigned)); /* array of min Hamming distances */

  swap_tree = complete_parse_nh(swap_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  for (i = 0; i < m; i++) {
    min_dist[i] = n; /* initialization to the nb of taxa */
  }
  min_num_moved = 4;
  e_index=6;

  /* calculation of the C and I matrices (see Brehelin/Gascuel/Martin) */
  update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
  update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);
  
  if(min_dist[e_index] != min_num_moved){
    fprintf(stderr,"TRANSFER Test 3 : Error : The min_dist of the internal branch is != %d (%d)\n",min_num_moved,min_dist[e_index]);
    exit(EXIT_FAILURE);
  }

  free_tree(swap_tree);

  swap_tree = complete_parse_nh(swap_tree2_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  for (i = 0; i < m; i++) {
    min_dist[i] = n; /* initialization to the nb of taxa */
  }
  min_num_moved = 5;
  e_index=6;

  /* calculation of the C and I matrices (see Brehelin/Gascuel/Martin) */
  update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
  update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);
  
  if(min_dist[e_index] != min_num_moved){
    fprintf(stderr,"TRANSFER Test 3 : Error : The min_dist of the internal branch is != %d (%d)\n",min_num_moved,min_dist[e_index]);
    exit(EXIT_FAILURE);
  }

  free_tree(swap_tree);


  for (i=0; i<m; i++) {
    free(c_matrix[i]);
    free(i_matrix[i]);
    free(hamming[i]);
  }
  free(c_matrix);
  free(i_matrix);
  free(hamming);
  free(min_dist);
  free_tree(ref_tree);

  fprintf(stderr,"TRANSFER Test 3 : OK\n");

  return EXIT_SUCCESS;
}


/**
   We test the TRANSFER Support for branches of the initial tree compared to another tree
*/
int test_transfer_4(){
  int seed = 103873987;
  int new_seed;
  /* int seed = 1038739; */
  Tree *ref_tree, *swap_tree;
  int i, i_edge;

  int n = 100;
  int n_swap = 10;
  int i_swap;
  double thresh = 0.05;
  int n_t = 10;
  int i_t;
  int collapsed_one, uncollapsed_terminal, collapsed_internal;
  prng_seed_bytes(&seed, sizeof(seed));

  /* We will perform the test for n_t simulated trees */
  for(i_t = 0; i_t < n_t; i_t++){
    new_seed = prng_get_int();
    prng_seed_bytes(&new_seed, sizeof(new_seed));
    ref_tree = gen_rand_tree(n, NULL);
    prng_seed_bytes(&new_seed, sizeof(new_seed));
    swap_tree = gen_rand_tree(n, NULL);
    
    /* Collapse branches */
    collapsed_internal = 0;
    do {
      collapsed_one = 0; /* flag that will be set to one as soon as we collapse one branch */
      uncollapsed_terminal = 0;
      for(i=0; i < ref_tree->nb_edges; i++) {
	if (ref_tree->a_edges[i]->brlen < thresh) {
	  if (ref_tree->a_edges[i]->right->nneigh == 1) { /* don't collapse terminal edges */
	    uncollapsed_terminal++;
	  }else{
	    collapse_branch(ref_tree->a_edges[i], ref_tree);
	    collapse_branch(swap_tree->a_edges[i], swap_tree);
	    collapsed_one = 1;
	    collapsed_internal++;
	    break; /* breaking the for so that we start again from the beginning because tree->a_edges has changed */
	  }
	}
      } /* end for */
    } while (collapsed_one);
    /* fprintf(stderr,"Collapsed %d branches\n",collapsed_internal); */
  
    int m = ref_tree->nb_edges;
    int max_branches_boot = ref_tree->nb_taxa*2-2;
    
    short unsigned** c_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of complements */
    short unsigned** i_matrix = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of cardinals of intersections */
    short unsigned** hamming = (short unsigned**) malloc(m*sizeof(short unsigned*)); /* matrix of Hamming distances */
    short unsigned* min_dist = (short unsigned*) malloc(m*sizeof(short unsigned)); /* array of min Hamming distances */
    
    for (i=0; i<m; i++) c_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
    for (i=0; i<m; i++) i_matrix[i] = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
    for (i=0; i<m; i++) hamming[i]  = (short unsigned*) malloc(max_branches_boot*sizeof(short unsigned));
    
    for(i_swap = 0; i_swap < n_swap; i_swap++){
      for (i = 0; i < m; i++) {
	min_dist[i] = n; /* initialization to the nb of taxa */
      }
      /* First we see if the min dist is 0 for all branches (it must be) */
      update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
      update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);
    
      for(i_edge=0; i_edge<ref_tree->nb_edges;i_edge++){
	if(min_dist[i_edge] != 0){
	  /* fprintf(stderr,"TRANSFER Test 4 : Error : The min_dist of the internal branch is != 0 (%d)\n",min_dist[i_edge]); */
	  return(EXIT_FAILURE);
	}
      }
      /* fprintf(stderr,"TRANSFER Test 4/1 : OK : The min_dist of the internal branch is == 0 (%d)\n",min_dist[i_edge]); */
    
      /* Then we exchange n_move taxa names from left to right of the edge */
      int edge = rand_to(ref_tree->nb_edges);
      int d = swap_tree->a_edges[edge]->topo_depth;
      int n_move = rand_to(d);
      /* fprintf(stderr,"\tWill swap %d taxa from left to right of branch %d (depth=%d)\n",n_move,edge,d); */
      int* left_taxa  = sample_taxa(swap_tree, n_move, swap_tree->a_edges[edge]->right, swap_tree->a_edges[edge]->left);
      int* right_taxa = sample_taxa(swap_tree, n_move, swap_tree->a_edges[edge]->left , swap_tree->a_edges[edge]->right);
      int i_move;
      for(i_move=0; i_move < n_move; i_move++){
	/* fprintf(stderr,"\tMoving %s <-> %s\n",swap_tree->a_nodes[left_taxa[i_move]]->name, swap_tree->a_nodes[right_taxa[i_move]]->name); */
	char *tmp;
	tmp = swap_tree->a_nodes[left_taxa[i_move]]->name;
	swap_tree->a_nodes[left_taxa[i_move]]->name = swap_tree->a_nodes[right_taxa[i_move]]->name;
	swap_tree->a_nodes[right_taxa[i_move]]->name = tmp;
      }
    
      for (i = 0; i < m; i++) {
	min_dist[i] = n; /* initialization to the nb of taxa */
      }
    
      update_all_i_c_post_order_ref_tree(ref_tree, swap_tree, i_matrix, c_matrix);
      update_all_i_c_post_order_boot_tree(ref_tree, swap_tree, i_matrix, c_matrix, hamming, min_dist);
    
      /* fprintf(stderr,"\tTRANSFER Test 4 : The min_dist of the internal branch is %d\n",min_dist[edge]); */
    
      if(min_dist[edge] > n_move*2 ){
	fprintf(stderr,"TRANSFER Test 4 : Error : The min_dist of the internal branch is > 2*%d (%d)\n",n_move,min_dist[edge]);
	return(EXIT_FAILURE);
      }

      /* We leave the tree as it was at the beginning */
      for(i_move=0; i_move < n_move; i_move++){
	char *tmp;
	tmp = swap_tree->a_nodes[left_taxa[i_move]]->name;
	swap_tree->a_nodes[left_taxa[i_move]]->name = swap_tree->a_nodes[right_taxa[i_move]]->name;
	swap_tree->a_nodes[right_taxa[i_move]]->name = tmp;
      }

      free(left_taxa);
      free(right_taxa);
    }
    
    for (i=0; i<m; i++) {
      free(c_matrix[i]);
      free(i_matrix[i]);
      free(hamming[i]);
    }
    free(c_matrix);
    free(i_matrix);
    free(hamming);
    free(min_dist);
    free_tree(ref_tree);
    free_tree(swap_tree);
  }
  return EXIT_SUCCESS;
}

int test_randomtree(){
  srand(time(NULL));
    char *ref_tree_string = "((a:1,b:1):1,e:1,(c:1,d:1):1);"; 
    char** taxname_lookup_table = NULL;
    Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
    int nbtrees = 10;
    int i,t;
    for(i=0;i<nbtrees;i++){
      Tree* rand_tree = gen_random_tree(ref_tree);
      /* We test if the lookup table is ok for all random tree*/
      for(t=0;t<rand_tree->nb_nodes;t++){
	Node *n = rand_tree->a_nodes[t];
	if(n->nneigh==1){
	  int ref_lookid = get_tax_id_from_tax_name(n->name,ref_tree->taxname_lookup_table, ref_tree->nb_taxa);
	  int rand_lookid= get_tax_id_from_tax_name(n->name,rand_tree->taxname_lookup_table, rand_tree->nb_taxa);
	  if(ref_lookid != rand_lookid){
	    fprintf(stderr,"Random tree test error: tax id in lookup table is : %d and should be : %d\n",rand_lookid,ref_lookid);
	    free_tree(rand_tree);
	    free_tree(ref_tree);
	    return EXIT_FAILURE;
	  }
	}
	/* For each Edge we will test the hashtables */
	int e;
	for(e=0;e<rand_tree->nb_edges;e++){
	  id_hash_table_t * h = rand_tree->a_edges[e]->hashtbl[1];
	  id_hash_table_t * h2 = create_id_hash_table(rand_tree->nb_taxa);
	  id_hash_table_t * h3 = create_id_hash_table(rand_tree->nb_taxa);
	  test_fill_hashtable_post_order(rand_tree->a_edges[e]->left,rand_tree->a_edges[e]->right, rand_tree, h2);
	  test_fill_hashtable_post_order(rand_tree->a_edges[e]->right,rand_tree->a_edges[e]->left, rand_tree, h3);
	  
	  if(!equal_id_hashtables(h,h2) && !equal_id_hashtables(h,h3)){
	  /* if(!equal_or_complement_id_hashtables(h,h2,rand_tree->nb_taxa)){ */
	    fprintf(stderr,"Random tree test error: hashtables are not consistent with the lookup table\n");
	    print_id_hashtable(stderr, h, rand_tree->nb_taxa);
	    print_id_hashtable(stderr, h2, rand_tree->nb_taxa);
	    print_id_hashtable(stderr, h3, rand_tree->nb_taxa);
	    
	    free_tree(rand_tree);
	    free_tree(ref_tree);
	    free_id_hashtable(h2);
	    return EXIT_FAILURE;	    
	  }
	  free_id_hashtable(h2);
	}
      }
      free_tree(rand_tree);
    }
    free_tree(ref_tree);

    fprintf(stderr,"Random tree Test: OK\n");

    return(EXIT_SUCCESS);
}

int test_id_hash_table_shuffle(){
  ntax = 1000;

  id_hash_table_t * h;
  id_hash_table_t * h2;
  int i=0;
  int total = 0;
  int total_expect = 0;
  h = create_id_hash_table(ntax);

  /* We will set the tax to 1 randomly */
  for(i = 0;i<ntax;i++){
    if(unif()<0.5){
      add_id(h,i);
      total_expect++;
    }
  }

  h2 = suffle_hash_table(h, ntax);

  if(equal_id_hashtables(h,h2)){
    fprintf(stderr,"Hashtable shuffle test error: the shuffled hash table is equal to the original one\n");
    free_id_hashtable(h);
    free_id_hashtable(h2);
    return EXIT_FAILURE;
  }

  for(i=0;i<ntax;i++){
    if(lookup_id(h2,i)){
      total++;
    }
  }
  if(total!=total_expect){
    fprintf(stderr,"Hashtable shuffle test error: the shuffled hash table has a different number of taxa than the original one: %d != %d\n",total,total_expect);
    free_id_hashtable(h);
    free_id_hashtable(h2);
    return EXIT_FAILURE;
  }
  free_id_hashtable(h2);
  free_id_hashtable(h);
  fprintf(stderr,"Hashtable shuffle Test: OK\n");
  return EXIT_SUCCESS;
}

int test_id_hash_table(){
  id_hash_table_t * h;
  id_hash_table_t * h2;
  id_hash_table_t * h3;
  /* Essai avec 125 taxons : Pour tester les 
     derniers bits à gauche: Il devrait y en avoir
     3 qui restent.
     2 chunks de 64 bits
   */
  id_hash_table_t * h4;
  id_hash_table_t * h5;

  ntax=5;
  h = create_id_hash_table(5);
  h2 = create_id_hash_table(5);
  h3 = create_id_hash_table(5);
  add_id(h,0);
  add_id(h,2);
  
  add_id(h2,1);
  add_id(h2,3);
  add_id(h2,4);

  add_id(h3,1);
  add_id(h3,3);
  add_id(h3,4);

  fprintf(stderr,"\t hashtable 1: ");
  print_id_hashtable(stderr, h, 5);
  fprintf(stderr,"\t hashtable 2: ");
  print_id_hashtable(stderr, h2, 5);
  fprintf(stderr,"\t hashtable 3: ");
  print_id_hashtable(stderr, h3, 5);

  if(equal_id_hashtables(h,h2)){
    fprintf(stderr,"Hash table Test error: the two hash tables must be different\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  if(!equal_or_complement_id_hashtables(h,h2,5)){
    fprintf(stderr,"Hash table Test error: the two hash tables should be equal or complement\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  if(!complement_id_hashtables(h,h2,5)){
    fprintf(stderr,"Hash table Test error: the two hash tables should be complement\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  if(!equal_id_hashtables(h2,h3)){
    fprintf(stderr,"Hash table Test error: the two hash tables should equal\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  if(!equal_or_complement_id_hashtables(h2,h3,5)){
    fprintf(stderr,"Hash table Test error: the two hash tables should equal\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  if(complement_id_hashtables(h2,h3,5)){
    fprintf(stderr,"Hash table Test error: the two hash tables should not be complement\n");
    print_id_hashtable(stderr, h, 5);
    print_id_hashtable(stderr, h2, 5);
    return EXIT_FAILURE;
  }

  ntax=125;
  h4 = create_id_hash_table(125);
  h5 = create_id_hash_table(125);
  int i=0;
  for(i=0;i<125;i+=2){
    add_id(h4,i);
    if(i<124)
      add_id(h5,i+1);
  }
  /* add_id(h4,124); */
  /* add_id(h5,124); */

  fprintf(stderr,"\t hashtable 4: ");
  print_id_hashtable(stderr, h4, 125);
  fprintf(stderr,"\t hashtable 5: ");
  print_id_hashtable(stderr, h5, 125);

  if(!equal_or_complement_id_hashtables(h4,h5,125)){
    fprintf(stderr,"Hash table Test error: the two hash tables should be equal or complement\n");
    print_id_hashtable(stderr, h4, 125);
    print_id_hashtable(stderr, h5, 125);
    return EXIT_FAILURE;
  }
  
  if(!complement_id_hashtables(h4,h5,125)){
    fprintf(stderr,"Hash table Test error: the two hash tables should be complement\n");
    print_id_hashtable(stderr, h4, 125);
    print_id_hashtable(stderr, h5, 125);
    return EXIT_FAILURE;
  }  

  if(equal_id_hashtables(h4,h5)){
    fprintf(stderr,"Hash table Test error: the two hash tables should not be equal\n");
    print_id_hashtable(stderr, h4, 125);
    print_id_hashtable(stderr, h5, 125);
    return EXIT_FAILURE;
  }

  fprintf(stderr,"Hash table Test: OK\n");
  return(EXIT_SUCCESS);
}

int test_stat_proba(){
  int expected=5000;
  int result = 0; 
  int i;
  for(i=0; i < 10000; i++){
    if(proba(0.5)){
      result++;
    }
  }
  /* We accept 5% error 250/5000 */
  if(abs(result-expected)>250){
    printf("Test proba : error - %d != %d\n",result,expected);
    return(EXIT_FAILURE);
  }
  printf("Test stat_proba : OK\n");
  return(EXIT_SUCCESS);
}

int test_qnorm(){
  int i;  
  double alphas[14] = {0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99};
  double expect[14] = {-2.32634787404084075746,
		      -2.05374891063182252182,
		      -1.64485362695147263601,
		      -1.28155156554460081253,
		      -0.84162123357291418468,
		      -0.52440051270804066696,
		      -0.25334710313579977825, 
		      0.00000000000000000000,
		      0.25334710313579977825,
		      0.52440051270804066696,
		      0.84162123357291440673,
		      1.28155156554460081253,
		      1.64485362695147152579,
		      2.32634787404084075746};
  double res;

  for(i=0; i <14; i++){
    res = qnorm(alphas[i], 0, 1);
    if(expect[i] != res){
      printf("Test qnorm : error - %f != %f\n",expect[i],res);
      return(EXIT_FAILURE);
    }
    /*printf("%f = %1.30f | %1.30f (%s)\n",alphas[i],expected[i],res,(expected[i]==res)?"true":"false");*/
  }
  printf("Test qnorm : OK\n");
  return(EXIT_SUCCESS);
}

int test_pnorm(){
  int i;  
  double q[5] = {-2,-1,0,1,2};
  double expect[5] = {
    0.022750131948179212055,
    0.15865525393145704647,
    0.5,
    0.84134474606854292578,
    0.97724986805182079141
  };

  double res;
  double relative_error;
  double accepted_error = 0.000000000000001;
  for(i=0; i <2; i++){
    res = pnorm(q[i]);
    relative_error = fabs((res - expect[i]));
    if(expect[i] != res && relative_error > accepted_error){
      printf("Test qnorm : error - %1.20f != %1.20f\n",expect[i],res);
      return(EXIT_FAILURE);
    }
  }
  printf("Test pnorm : OK\n");
  return(EXIT_SUCCESS);
}

int test_unif(){
  int seed = 25684;
  double expected[7] = {0.614942,0.295840,0.981761,0.359667,0.436287,0.827348,0.813658};
  double result;
  int i;

  prng_seed_bytes(&seed, sizeof(seed));

  for(i = 0; i < 7; i++){
    result = unif();
    if((int)round(expected[i]*1000000) != (int)round(result*1000000)){
      printf("Test unif : error - %d != %d\n",(int)(expected[i]*1000000),(int)(result*1000000));
      return(EXIT_FAILURE);
    }
  }
  printf("Test unif : OK\n");
  return(EXIT_SUCCESS);
}

#define TEST_MAX_INT 124
int test_rand_to(){
  int nb_simu = 1000000;
  double expected_nb = nb_simu/TEST_MAX_INT;
  double threshold = 0.1;
  double res_nb [TEST_MAX_INT];
  int i;
  int result;

  prng_seed_time();

  for(i = 0; i < TEST_MAX_INT; i++){
    res_nb[i] = 0;
  }

  for(i = 0; i < nb_simu; i++){
    result = rand_to(TEST_MAX_INT);
    if(result >= TEST_MAX_INT){
      printf("Test rand_to : integer %d > %d\n",result,TEST_MAX_INT-1);
    }
    res_nb[result]++;
  }
  /* Test for frequency of each int */
  for(i=0;i<TEST_MAX_INT;i++){
    if(fabs(res_nb[i]-expected_nb) > threshold*expected_nb){
      printf("Test rand_to : frequency error - freq %d = %f != %f \n",i,res_nb[i],expected_nb);
      return(EXIT_FAILURE);
    }
  }

  printf("Test rand_to: OK (~1/TEST_MAX_INT of each nt)\n");
  return(EXIT_SUCCESS);
}

int test_sum(){
  double array[10] = {1,2,3,4,5,6,7,8,9,10};
  double result = sum(array,10);
  double exp = 55;
  if(result!=exp){
    printf("Test sum : error - Sum %f != %f\n",result,exp);
  }
  printf("Test sum: OK\n");
  return(EXIT_SUCCESS);  
}

int comp_int(const void * elem1, const void * elem2){
  int f = *((int*)elem1);
  int s = *((int*)elem2);
  if (f > s) return  1;
  if (f < s) return -1;
  return 0;
}

/* Test sampling function */
int test_sample(){
  int length = 10000;
  int nbsamp = 500;
  int * array = malloc(length * sizeof(int));
  int * sampled;
  int * sampled2;
  int found, found2;
  int i = 0, j = 0;
  int duplicate = 0;

  for(i=0;i<length;i++){
    array[i] = i;
  }
  sampled  = sample(array, length, nbsamp, 0);
  sampled2 = sample(array, length, nbsamp, 1);

  /* Test if all the output values are in the original array */
  for(i = 0; i < nbsamp; i++){
    found = 0;
    found2 = 0;
    for(j = 0; j < length; j++){
      if(sampled[i] == array[j]){
	found = 1;
      }
      if(sampled2[i] == array[j]){
	found2 = 1;
      }
    }
    if(!found){
      fprintf(stderr,"Test sample : error - The sampled value %d is not in the original data\n",sampled[i]);
      free(array);
      free(sampled);
      free(sampled2);
      return(EXIT_FAILURE);
    }
    if(!found2){
      fprintf(stderr,"Test sample : error - The sampled2 value %d is not in the original data\n",sampled2[i]);
      free(array);
      free(sampled);
      free(sampled2);
      return(EXIT_FAILURE);
    }
  }

  /* Test if the sampled arrays are different from the al arrays */
  int all_equals=1, all_equals2 = 1;
  for(i=1;i<nbsamp;i++){
    if(sampled[i] != array[i]){
      all_equals = 0;
    }
    if(sampled2[i] != array[i]){
      all_equals2 = 0;
    }
  }
  if(all_equals){
    fprintf(stderr,"Test sample : error - The %dth sampled values are the same than the original values\n",nbsamp);
    free(array);
    free(sampled);
    free(sampled2);
    return(EXIT_FAILURE);
  }
  if(all_equals2){
    fprintf(stderr,"Test sample : error - The %dth sampled2 values are the same than the original values\n",nbsamp);
    free(array);
    free(sampled);
    free(sampled2);
    return(EXIT_FAILURE);
  }

  /* Test if no duplicate values in the sample without replacement */
  qsort(sampled, nbsamp, sizeof(int), comp_int);
  for(i=1;i<nbsamp;i++){
    if(sampled[i-1] == sampled[i]){
      fprintf(stderr,"Test sample: error - Duplicate values after sample: %d\n",sampled[i]);
      free(array);
      free(sampled);
      free(sampled2);
      return(EXIT_FAILURE);
    }
  }
  /* Test if duplicate values in the sample with replacement */
  qsort(sampled2, nbsamp, sizeof(int), comp_int);
  duplicate=0;
  for(i=1;i<nbsamp;i++){
    if(sampled2[i-1] == sampled2[i]){
      duplicate=1;
    }
  }
  if(!duplicate){
    fprintf(stderr,"Test sample: error - No duplicate values after sample with replacement\n");
    free(array);
    free(sampled);
    free(sampled2);
    return(EXIT_FAILURE);
  }

  free(array);
  free(sampled);
  free(sampled2);
  fprintf(stderr,"Test sample: OK\n");
  return(EXIT_SUCCESS);
}


int test_sample_from_counts(){
  int i;
  int input1[4] = {1,0,0,0};
  int input2[4] = {0,1,1,1};
  int input3[4] = {1,1,1,1};
  int input4[4] = {0,0,0,0};
  int input5[4] = {0,0,1,0};

  int* output1 = sample_from_counts(input1, 4, 1, 0);
  int* output2 = sample_from_counts(input2, 4, 3, 0);
  int* output3 = sample_from_counts(input3, 4, 2, 0);
  int* output4 = sample_from_counts(input4, 4, 1, 0);
  int* output5 = sample_from_counts(input5, 4, 2, 0);

  int sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0,sum5 = 0;

  for(i=0;i<4;i++){
    sum1 += output1[i];
    sum2 += output2[i];
    sum3 += output3[i];
    sum4 += output4[i];
    sum5 += output5[i];
    
    if(output4[i]!=0 || output5[i]!=0){
      fprintf(stderr,"Test sample from counts: error - The array must be 0 filled\n");
      free(output1);
      free(output2);
      free(output3);
      free(output4);
      free(output5);
      return(EXIT_FAILURE);
    }
  }

  if(sum1!=1 || sum2!=3 || sum3!=2 || sum4 != 0 || sum5 != 0){
      fprintf(stderr,"Test sample from counts: error - The counts do not sum to the expected total\n");
      free(output1);
      free(output2);
      free(output3);
      free(output4);
      free(output5);
      return(EXIT_FAILURE);
  }

  fprintf(stderr,"Test sample from counts: OK\n");
  return(EXIT_SUCCESS);
}


/**
   test of original tree:
     a   e   d      a   b   d         c   b   d         c     d               a  d
      \  |  /        \  |  /           \  |  /           \   /                 \/ 
       .-.-.    vs.   .-.-.   and vs.   .-.-.   and vs.   .--.-a   and vs.   c--.--b 
      /     \	     /     \ 	       /     \ 	          /   \                 | 
     b       c      e       c         e       a	         e     d                e
 */
int test_remove_taxon(){
  char *ref_tree_string = "((a:10,b:1.5)3:0.2,e:1,(c:1,d:1)2:1);"; 
  char *boot3_tree_string = "((c:1,e:1):1,a:1,b:1,d:1);";
  char *boot4_tree_string = "(c:1,a:1,b:1,d:1,e:1);";

  int i = 0;
  double sum_brlen = 0.0;

  char** taxname_lookup_table = NULL;
  Tree* ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  int tax_id = get_tax_id_from_tax_name("a", taxname_lookup_table, ref_tree->nb_taxa);
  sum_brlen=0.0;
  write_nh_tree(ref_tree,stdout);
  remove_taxon(tax_id,ref_tree);

  for(i=0;i<ref_tree->nb_edges;i++){
    sum_brlen += (ref_tree->a_edges[i]->brlen);
  }
  if(sum_brlen != 5.7){
    fprintf(stderr,"Test remove taxon: error - The sum of br len is %f, and should be 5.7\n",sum_brlen);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }

  if(ref_tree->nb_nodes != 6){
    fprintf(stderr,"Test remove taxon: error - The number of nodes is %d, and should be 6\n",ref_tree->nb_nodes);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }
  if(ref_tree->nb_taxa != 4){
    fprintf(stderr,"Test remove taxon: error - The number of taxa is %d, and should be 4\n",ref_tree->nb_taxa);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }

  for(i=0;i<ref_tree->nb_nodes;i++){
    if(ref_tree->a_nodes[i]->nneigh==1 && strcmp(ref_tree->a_nodes[i]->name,"a") == 0){
      fprintf(stderr,"Test remove taxon: error - The original taxon \"a\" is still present after its removal\n");
      free_tree(ref_tree);
      free(taxname_lookup_table);
      return(EXIT_FAILURE);
    }
  }
  write_nh_tree(ref_tree,stdout);
  free(taxname_lookup_table);
  taxname_lookup_table=NULL;
  free_tree(ref_tree);

  /* On essaie avec un arbre multifurcation */
  ref_tree = complete_parse_nh(boot4_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  tax_id = get_tax_id_from_tax_name("a", taxname_lookup_table, ref_tree->nb_taxa);
  remove_taxon(tax_id,ref_tree);
  if(ref_tree->nb_nodes != 5){
    fprintf(stderr,"Test remove taxon: error - The number of nodes is %d, and should be 5\n",ref_tree->nb_nodes);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }
  if(ref_tree->nb_taxa != 4){
    fprintf(stderr,"Test remove taxon: error - The number of taxa is %d, and should be 4\n",ref_tree->nb_taxa);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }

  for(i=0;i<ref_tree->nb_nodes;i++){
    if(ref_tree->a_nodes[i]->nneigh==1 && strcmp(ref_tree->a_nodes[i]->name,"a") == 0){
      fprintf(stderr,"Test remove taxon: error - The original taxon \"a\" is still present after its removal\n");
      free_tree(ref_tree);
      free(taxname_lookup_table);
      return(EXIT_FAILURE);
    }
  }
  write_nh_tree(ref_tree,stdout);
  free(taxname_lookup_table);
  taxname_lookup_table=NULL;
  free_tree(ref_tree);

  /* On essaie avec un autre arbre multifurcation */
  ref_tree = complete_parse_nh(boot3_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  tax_id = get_tax_id_from_tax_name("a", taxname_lookup_table, ref_tree->nb_taxa);
  remove_taxon(tax_id,ref_tree);
  if(ref_tree->nb_nodes != 6){
    fprintf(stderr,"Test remove taxon: error - The number of nodes is %d, and should be 6\n",ref_tree->nb_nodes);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }
  if(ref_tree->nb_taxa != 4){
    fprintf(stderr,"Test remove taxon: error - The number of taxa is %d, and should be 4\n",ref_tree->nb_taxa);
    free_tree(ref_tree);
    free(taxname_lookup_table);
    return(EXIT_FAILURE);
  }

  for(i=0;i<ref_tree->nb_nodes;i++){
    if(ref_tree->a_nodes[i]->nneigh==1 && strcmp(ref_tree->a_nodes[i]->name,"a") == 0){
      fprintf(stderr,"Test remove taxon: error - The original taxon \"a\" is still present after its removal\n");
      free_tree(ref_tree);
      free(taxname_lookup_table);
      return(EXIT_FAILURE);
    }
  }

  write_nh_tree(ref_tree,stdout);
  free(taxname_lookup_table);
  taxname_lookup_table=NULL;
  free_tree(ref_tree);
  

  ref_tree = complete_parse_nh(ref_tree_string, &taxname_lookup_table); /* sets taxname_lookup_table en passant */
  tax_id = get_tax_id_from_tax_name("e", taxname_lookup_table, ref_tree->nb_taxa);
  remove_taxon(tax_id,ref_tree);
  write_nh_tree(ref_tree,stdout);

  free_tree(ref_tree);
  free(taxname_lookup_table); /* which is a (char**) */
  fprintf(stderr,"Test remove taxon: OK\n");
  return(EXIT_SUCCESS);

}

int test_hashmap(){

  int j,k,div,mod;
  int total=20000;
  char* array[total];
  int min_char=65;
  int max_char=90;
  int base=max_char-min_char+1;
  /* Wi fill the array of string with strings AAAAAAAAA then BAAAAAAAA, etc.*/
  for(j=0;j<total;j++){
    array[j] = malloc(10*sizeof(char));
    div=j;
    mod=0;
    for(k=0;k<9;k++){
      if(div==0){
	array[j][k]=65;
      }else{
	mod=div%base;
	div=div/base;
	array[j][k] = mod+65;
      }
    }
    array[j][9] = '\0';
  }

  map_t h = hashmap_new();

  int i;
  for(i=0;i<total;i++){
    int *val = malloc(sizeof(int));
    *val=i;
    hashmap_put(h, array[i], val);
  }

  for(i=0;i<total;i++){
    int *val;
    hashmap_get(h, array[i],(void**) &val);
    if(*val!=i){
      fprintf(stderr,"Test hashmap: error - Key [%s] The Value from the HashMap (%d) is different from the original one (%d)\n",array[i],*val,i);
      return(EXIT_FAILURE);
    }
  }
  free_taxid_hashmap(h);
  
  fprintf(stderr,"Test hashmap: OK\n");  
  return(EXIT_SUCCESS);
}

int main(int arbc, char** argv){
  srand(time(NULL)); /* seeding the random generator */
  
  int exit_code = test_id_hash_table();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_id_hash_table_shuffle();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_swap_branches();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_randomtree();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_classical_bootstrap();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_stat_proba();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_qnorm();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_pnorm();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_unif();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_rand_to();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_sum();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_sample();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_sample_from_counts();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_remove_taxon();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_hashmap();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_transfer_1();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_transfer_2();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_transfer_3();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }

  exit_code = test_transfer_4();
  if(exit_code != EXIT_SUCCESS){
    return(exit_code);
  }
  
  return(exit_code);
}

