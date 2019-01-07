#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

float tab_pcontact[1024][1024];
float pcontact;
int iteration=1;
int nbre_de_coupe=0;
int max_i1=0;
int max_i2=0;
int max_j1=0;
int max_j2=0;
int start;
int end;
int max_start=0;
int max_end=0;
int idx;

idx=0;

#define MONOMER '_'
#define MAX_ITERATION 32

float max_coeff_matthews=0.;
int size_pu;
int min_size_pu;

int LIMITSIZESS2;
int LIMITSIZEPU;
int MAXSIZEPU;

int MAXR2;

int MIN_SIZE_PU;
int MAX_SIZE_PU;

float D0;
float DELTA;
int ONLYSS2;

int id_pu=0;

int nbre_pu=0;
int new_nbre_pu=0;
int pu[32][128][2];

int best_pu=-1;
int current_pu;
float CI;


int tab_decoupe[1024];
int tab_ss2[1024];

/* Pruning and Homogeneity variable */
int PRUNING;
float CUTOFF_PRUNING;
double H_PU1;
double H_PU2;
double H_PU3;
int start_pu;
int end_pu;

int main (int argc, char *argv[])
{
	/* Protype des fonctions */
	void parse_pdb(char *, char *,char *); /* Protype de la fonction parse_pdb */
	void parse_dssp(char *); /* Protoype de la fonction parse_dssp */
	void simple_cutting(void);
	void double_cutting(void);
	double homogeneity(int ,int);
	void save_pu(void);
	void mutual_information(void);
	void test(void);

	char NAME_PDB_FILE[1024];
	char NAME_DSSP_FILE[1024];
	char NAME_LIMITSIZESS2[5];
	char NAME_LIMITSIZEPU[5];
	char NAME_MAXSIZEPU[5];
	char NAME_R2[5];
	char NAME_D0[5];
	char NAME_DELTA[5];
	char NAME_ONLYSS2[5];
	char NAME_PRUNING[5];
	char NAME_CUTOFF_PRUNING[5];

	char code_pdb[5];
	char *p_name_pdb_file;
	char *p_name_dssp_file;
	char *p_code_pdb;

	int continuation;

	int x,y,z;

	int i;

	/*Recuperation des arguments */
	strcpy (NAME_PDB_FILE,argv[1]);      //argument 1: PDB file
	strcpy (NAME_DSSP_FILE,argv[2]);     //argument 2: DSSP file
	strcpy (NAME_R2,argv[3]);            //argument 3: R2 max
	strcpy (NAME_LIMITSIZESS2,argv[4]);  //argument 4: limit size of secondary structure wo cutting
	strcpy (NAME_LIMITSIZEPU ,argv[5]);  //argument 5: inferior limit size of PU
	strcpy (NAME_MAXSIZEPU ,argv[6]);    //argument 6: if PU have a size > limit, peeling continu
	strcpy (NAME_D0, argv[7]);           //argument 7: D0 for logistic function
	strcpy (NAME_DELTA, argv[8]);        //argument 8: DELTA for logistic function
	strcpy (NAME_ONLYSS2, argv[9]);      //argument 9: Only Calpha include in SS2 for computation of PU (other Calpha have their probability divide by 10
	strcpy (NAME_PRUNING, argv[10]);    //argument 10: Pruning tree ? (yes|no -> 1|0 )
	strcpy (NAME_CUTOFF_PRUNING, argv[11]);   //argument 11: cutoff value of CI (HI)

	/*Conversion char en entier */
	sscanf (NAME_R2,          "%d" ,&MAXR2);
	sscanf (NAME_LIMITSIZESS2,"%d" ,&LIMITSIZESS2);
	sscanf (NAME_LIMITSIZEPU, "%d" ,&LIMITSIZEPU);
	sscanf (NAME_MAXSIZEPU   ,"%d" ,&MAXSIZEPU);
	sscanf (NAME_PRUNING    ,"%d", &PRUNING);
	sscanf (NAME_ONLYSS2,"%d" ,&ONLYSS2);

	/* Convertion char in float */
	sscanf (NAME_D0,"%f" ,&D0);
	sscanf (NAME_DELTA,"%f" ,&DELTA);
	sscanf (NAME_CUTOFF_PRUNING, "%f", &CUTOFF_PRUNING);

	MIN_SIZE_PU=LIMITSIZEPU;
	MAX_SIZE_PU=MAXSIZEPU;
	//printf("args: MAXR2:%d LIMITSIZE:%d MAXSIZEPU:%d NAME_D0:%f NAME_DELTA:%f NAME_ONLYSS2:%d\n",MAXR2,LIMITSIZESS2,MAX_SIZE_PU,D0,DELTA,ONLYSS2);
	//printf ("MAXR2:%d LIMITSIZE:%d\n",MAXR2,LIMITSIZESS2);
	//printf ("MAXSIZEPU:%d \n",MAX_SIZE_PU);
	/* Affichage des arguments*/
	//strncpy(code_pdb,NAME_PDB_FILE,5);
	//code_pdb[5]='\0';
	p_name_pdb_file=NAME_PDB_FILE;
	p_name_dssp_file=NAME_DSSP_FILE;
	p_code_pdb=code_pdb;

	/* Mise a zero du tableau decoupee */
	for (i=start ; i <= end ; i++)
	  {
		tab_decoupe[i]=1;
	  }

	/* Parsepdb */
	parse_pdb(p_code_pdb,p_name_pdb_file, p_name_dssp_file);
	/* Parsedssp */
	parse_dssp(p_name_dssp_file);

	nbre_pu=0; /* A la premiere iteration il y a une seule pu*/
	iteration=1;



	while (iteration < 32)
	  {
		/* Affectation */
		nbre_pu=new_nbre_pu;
		/* Initialisation */
		new_nbre_pu=0;
		best_pu=-1;
		for (x=0 ; x <= nbre_pu ; x++)
		  {
			current_pu=x;
			/*Assignation du nbre de pu */
			if (iteration > 0)
			  {
				start=pu[iteration-1][x][0];
				end=pu[iteration-1][x][1];
			  }
			/* Verification de la taille de la pu */
			size_pu=end-start;
			min_size_pu=MIN_SIZE_PU;
			if (size_pu >= min_size_pu)
			  {
				/*Simple decoupe */
				simple_cutting();
				/* double decoupe*/
				double_cutting();
			  }
		  }
		if (best_pu == -1) {exit(1);}


		/*Sauvegarde des unites coupees*/
		save_pu();

		/* Verification de la taille max des pu si variable minmax_pu est differente de 0*/
		continuation=0;
		if (MAX_SIZE_PU != 0 )
		{
			continuation=0;
			for (x=0 ; x <= nbre_pu ; x++)
			{
				current_pu=x;
				/*Assignation du nbre de pu */
				if (iteration > 0)
				{
					start=pu[iteration-1][x][0];
					end=pu[iteration-1][x][1];
				}
				/* Verification de la taille de la pu */
				size_pu=end-start;
				if (size_pu > MAX_SIZE_PU)
				{
					continuation=1;
				}
			}
			if (continuation == 0 ) {exit(1);}
		}

		/* Verification de la valeur du CI/HI (homogenity idx) des PUs crees et elegage le cas echeant */
		if (PRUNING != 0)
		{
			continuation=0;
			if (nbre_de_coupe==1)
			{
				/* Premiere PU */
				H_PU1=homogeneity(max_start,max_i1);
				/* Deuxieme PU */
				H_PU2=homogeneity(max_i2,max_end);
				/* Test valeur du CI/HI*/
				if (H_PU1 >= CUTOFF_PRUNING || H_PU2 >= CUTOFF_PRUNING )
				{
					//printf("Continuation 1 %f  >= %f ou  %f  >= %f !\n",H_PU1, CUTOFF_PRUNING,H_PU2,CUTOFF_PRUNING);
					continuation=1;
				}
			}
			else if (nbre_de_coupe==2)
			{
				/* Premiere PU */
				H_PU1=homogeneity(max_start,max_i1);
				/* Deuxieme PU */
				H_PU2=homogeneity(max_i2,max_j1);
				/* Toisieme PU */
				H_PU3=homogeneity(max_j2,max_end);
				/* Test valeur du CI/HI*/
				if (H_PU1 >= CUTOFF_PRUNING || H_PU2 >= CUTOFF_PRUNING || H_PU3 >= CUTOFF_PRUNING)
				{
					//printf("Continuation 2 %f >= %f ou  %f >= %f  %f >= %f !\n",H_PU1, CUTOFF_PRUNING,H_PU2,CUTOFF_PRUNING,H_PU3,CUTOFF_PRUNING);
					continuation=1;
				}
			}
			if (continuation == 0 ) {exit(1);}
		}



		/*mutual_information();*/
		mutual_information();
		max_coeff_matthews=0.;
		iteration++;

	  }



	//printf ("END PU PARSING\n");
}
/*
   1--------10-------20--------30--------40--------50--------60
   +--------+--------+---------+---------+---------+---------+
   ATOM      2  CA  ALA     1      11.846  49.175  18.125  1.00 23.06       1SG   3
 */

/****************************************************************************/
/*                                                                          */
/*  FONCTION OUVERTURE ET TRAITEMENT DU PDB ET CALCUL MATRICE DE CONTACT   */
/*                                                                          */
/****************************************************************************/
void parse_pdb (char *p_code_pdb, char *p_name_pdb_file, char *p_name_dssp_file)
{

	/* Variables lecture fichier */
	char line[85];
	char line_out[80];
	FILE *pdb_file;
	FILE *dssp_file;
	char chain='x';
	char at, ch, aa3[4], path[128], str[128], buf[128], atname[4],aaname[3];
	float x=0;
	float y=0;
	float z=0;
	int numat;
	int c=1;/* variable contenant le nombre de chaine */
	int tab_x[1024],tab_y[1024],tab_z[1024];
	char tab_chain[1024];

	/* Variable calcul des distances */
	int i;
	int j;

	int dx2,dy2,dz2;
	float p;
	double dt,dt2;
	float tmp;

	/* Varaibles calculs des probabilite de contacts */
	float mean_pcontact;
	float co;
	float pcontactaa;

	float pccoaa,plcoaa;

	int tab_contact[1024][1024];

	int delta_ij;
	int nb_p;

	idx = 0;

	int n;
	int L;
	int L2;

	/* Varables decriture de fichier */
	char *file_dist_contact_mtx ="file_dist_contact.mtx";
	char *file_proba_contact_mtx   ="file_proba_contact.mtx";
	char *file_dist_contact_mat ="file_dist_contact.mat";
	char *file_proba_contact_mat   ="file_proba_contact.mat";

	char *file_position_contactorder="file_position_contactorder";
	char *file_position_inv_contactorder="file_position_inv_contactorder";
	char *file_position_dist_contact="file_position_dist_contact";
	char *file_position_proba_contact="file_position_proba_contact";


	char *file_ca_coo ="file_ca_coo.pdb";

	//printf("%s %s %s LIMITSIZE:%d\n",p_code_pdb,p_name_pdb_file, p_name_dssp_file, LIMITSIZESS2);
	/*Ouverture du fichier file_ca_coo*/
	FILE* pfile_ca_coo = fopen(file_ca_coo, "w");
	pdb_file=fopen(p_name_pdb_file, "r"); /*Ouvrir en lecture*/

	if (pdb_file == NULL)
	{
		printf ("Impossible d'ouvir le fichier PDB !! %s\n", p_name_pdb_file);
		exit(1);
	}

	while (!feof(pdb_file))
	{
		fgets(line,85,pdb_file);
		if (strncmp("ATOM", line, 4)) continue; /* skip non ATOM */

		/* Execute */
		at = line[13];

		/* Chain extraction */
		char useless_str[6];
		sscanf(line, "%6s %6d %4s", useless_str, &numat, atname);

		//sscanf(&line[6] ,"%d", &numat );
		//sscanf(&line[13],"%3s", atname);

		if (strncmp("CA",atname,2)) continue; /* skip non CA */
		sscanf(&line[15],"%3s", aaname);

		/* Si c'est la première fois que l'on rencontre le champs chaine */
		if (chain == 'x')
		{
			if (line[21]==' ')
			{
				ch = MONOMER;
			}
			else
			{
				ch = line[21];
				if (chain==MONOMER) chain=ch; /* if chain=='_', keep first chain */
			}
		}

		//printf("%i %s", idx,line);
		/* Extraction coordonnes Atom Calpha */
		line[54]=' '; sscanf(&line[46], "%f", &z);
		line[46]=' '; sscanf(&line[38], "%f", &y);
		line[38]=' '; sscanf(&line[30], "%f", &x);

		/*printf("%d ATOM %d %-2s %3s ",idx,numat,atname,aaname);*/
		/*printf("%d %s %s %f %f %f \n",idx,aaname,atname,x,y,z);*/
		tab_chain[idx]=ch;
		tab_x[idx]=x;
		tab_y[idx]=y;
		tab_z[idx]=z;
		idx++; /* Ajouter 1 compteur idx */
		strncpy(line_out,line,54);
		fprintf(pfile_ca_coo,"%s\n",line_out);
	}
	idx--; // elimination de la derniere incremetation car elle est vide
	fclose(pdb_file); /*Fermer le fichier */
	fclose(pfile_ca_coo);
	/* Calcul distance intercalpha dans un tableau*/
	pcontact=0;
	pccoaa=0.;
	plcoaa=0.;
	nb_p=0;

	//printf("idx end:%i\n",idx);
	L=idx+1;
	L2=L*L;

	/*Ouverture des fichiers*/
	//FILE* pfile_dist_contact_mtx = fopen(file_dist_contact_mtx, "w");
	//FILE* pfile_proba_contact_mtx   = fopen(file_proba_contact_mtx,   "w");
	//FILE* pfile_dist_contact_mat = fopen(file_dist_contact_mat, "w");
	FILE* pfile_proba_contact_mat = fopen(file_proba_contact_mat, "w");


	//FILE* pfile_position_contactorder=fopen(file_position_contactorder,"w");
	//FILE* pfile_position_inv_contactorder=fopen(file_position_inv_contactorder,"w");
	//FILE* pfile_position_dist_contact=fopen(file_position_dist_contact,"w");
	//FILE* pfile_position_proba_contact=fopen(file_position_proba_contact,"w");

	//fprintf(pfile_dist_contact_mtx,"%d %d %d\n", L, L, L2);
	//fprintf(pfile_proba_contact_mtx,"%d %d %d\n", L, L, L2);

	for (i=0 ; i <= idx ; i++)
	{
		for (j=0 ; j <= idx ; j++)
		{

			/* Calcul de la distance inter atome*/
			dx2=pow((tab_x[i]-tab_x[j]),2);
			dy2=pow((tab_y[i]-tab_y[j]),2);
			dz2=pow((tab_z[i]-tab_z[j]),2);
			dt2=dx2+dy2+dz2;
			dt = sqrt (dt2);

			/*Ecriture dans les fichiers de sortie */
			//fprintf (pfile_dist_contact_mat,"%5.2f ",dt);
			//fprintf (pfile_dist_contact_mtx,"%3d %3d %5.2f\n",i+1,j+1,dt);

			/*Calcul des proba selon un modele logistique*/
			tmp=dt-D0;
			tmp=tmp/DELTA;
			p=1/(1+exp(tmp));
			/*Ecriture dans les fichiers de sortie*/
			fprintf(pfile_proba_contact_mat,"%5.3f ",p);
			//fprintf (pfile_proba_contact_mtx,"%3d %3d %5.3f\n",i+1,j+1,p);

			/*Calcul du p contact total*/
			pcontact=pcontact+p ;
			/* Calcul du pseudo contact order */
			tab_contact[i][j]=0;
			if (j != i)
			{
				if (dt < D0)
				{
					tab_contact[i][j]=1;
				}
				pccoaa= pccoaa+(p/abs(i-j));
				plcoaa=plcoaa+(p*abs(i-j));
				nb_p++;
			}
			//fprintf (pfile_position_dist_contact,"%d ",tab_contact[i][j]);
			/*Insertion dans tableau pour calcul*/
			tab_pcontact[i][j]=p;
			/* Calcul pcontact par aa */
			pcontactaa=pcontactaa+p;
		}
		//fprintf (pfile_dist_contact_mat,"\n ");
		fprintf (pfile_proba_contact_mat,"\n ");

		//fprintf (pfile_position_dist_contact,"\n");

		/* Ecriture */
		/* Ecriture du pccoaa et plcoaa*/
		/*Rinitialisation de pcontactaa*/
		pcontactaa=0;

		plcoaa=plcoaa/nb_p++;
		/* Imprimer dans fichier valeur de pccoaa et pcloaa */
		//fprintf (pfile_position_inv_contactorder,"%f\n",plcoaa);
		//fprintf (pfile_position_contactorder,"%f\n",pccoaa);
		/*Reenitialisation de pcco et plco*/
		pccoaa=0.;
		plcoaa=0.;
		nb_p=0;
	}

	/*Fermeture des fichiers de sortie*/
	//fclose (pfile_dist_contact_mat);
	fclose(pfile_proba_contact_mat);
	//fclose (pfile_dist_contact_mtx);
	//fclose (pfile_proba_contact_mtx);

	//fclose (pfile_position_inv_contactorder);
	//fclose (pfile_position_contactorder);
	//fclose (pfile_position_dist_contact);

	start=0;
	end=idx;
	pu[0][0][0]=start;
	pu[0][0][1]=end;
	iteration++;
	best_pu=0;
}

/****************************************************************************/
/*                                                                          */
/*  FONCTION OUVERTURE ET TRAITEMENT DU FICHIER DSSP                        */
/*                                                                          */
/****************************************************************************/
/*   0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
 *   0....;....10...;....20...;....30...;....40...;....50...;....60...;....70...;....80...;....90...;....100..;....110..;....120..;....130
 *   #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
 *
 * */
void parse_dssp (char *p_name_dssp_file)
{
	/*Variables */
	FILE *dssp_file;

	char line[200];
	char endline[1];
	int numaa,numaa2;
	char aaname[1];
	char ss2[2];
	int idx2=0;
	int i,j,k;

	int start1=0;
	int sizess2;
	//int LIMITSIZESS2 = 8;
	int startss2,endss2;
	int tab_u_ss2[128];

	for (i=0 ; i <= 199 ; i++){line[i]='X'; }
	for (i=0 ; i <= 1024; i++){tab_ss2[i]=0;}

	//printf ("LIMITSIZESS2:%d\n",LIMITSIZESS2);

	dssp_file=fopen(p_name_dssp_file, "r"); /*Ouvrir en lecture*/
	if (dssp_file == NULL) {printf ("Impossible d'ouvir le fichier %s\n", p_name_dssp_file);exit (1);}
	while (!feof(dssp_file))
	{
		fgets(line,200,dssp_file);
		//if (idx2 == 0) {idx2++; continue;} /* Skip first line */
		line[127]='\0';
		endline[0]=line[126];
		//printf ("CMP:%s %c\n",endline, line[126]);
		endline[1]='\0';
		if (isspace(endline[0])) {continue;} /* retourne un caractere different de 0 si le caracatere est un espace et skippe*/
		if (strcmp("-",endline)==0 ) {continue;}
		/* Execute */
		ss2[0]=line[16];
		ss2[1]='\0';
		//*ss2='\0';
		/* Chain extraction */

		sscanf(&line[4] ,"%d", &numaa );
		sscanf(&line[9] ,"%d", &numaa2);
		sscanf(&line[13],"%1s", aaname);

		//printf("%-5d |%s|fion \n", idx2, ss2);
		if (isspace(ss2[0]) == 0) {tab_ss2[idx2]=0;}
		if (strcmp("H",ss2) == 0) {tab_ss2[idx2]=1;} /* */
		if (strcmp("B",ss2) == 0) {tab_ss2[idx2]=2;}
		if (strcmp("E",ss2) == 0) {tab_ss2[idx2]=2;}
		if (strcmp("G",ss2) == 0) {tab_ss2[idx2]=1;}

		idx2++; /* Ajouter 1 compteur idx */
		ss2[0]=' ';
	}
	idx2--; //Elimination de la dernire incrementation qui est vide
	//printf("idx2:%i\n",idx2);
	fclose(dssp_file); /*Fermer le fichier */

	j=0;
	for (i=0 ; i <= idx2 ; i++)
	{
		//printf("%-5d |%d|\n",i,tab_ss2[i]);
		/* Si on a une ss2 et que start1 est == 0 */
		if (tab_ss2[i] != 0  && start1 == 0)
		{
			start1 = tab_ss2[i];
			startss2=i;
		}
		/* Si on change de type de structure secondaire */
		if (tab_ss2[i] != start1 && tab_ss2[i] != 0)
		{
			start1 = tab_ss2[i];
			endss2=i-1;
			sizess2=endss2-startss2;
			if (sizess2 <= LIMITSIZESS2 && sizess2 != 0)
			{
				tab_u_ss2[j]=startss2;
				j++;
				tab_u_ss2[j]=endss2;
				j++;
			}

		}
		/* Si apres une serie de structure secondaire on retrouve une position sans structure secondaire */
		if (tab_ss2[i] == 0 && start1 != 0)
		{
			start1=0;
			endss2=i-1;
			sizess2=endss2-startss2;
			if (sizess2 <= LIMITSIZESS2 && sizess2 != 0)
			{
				tab_u_ss2[j]=startss2;
				j++;
				tab_u_ss2[j]=endss2;
				j++;
			}

		}
	}
	/* Initialisation du tableau autorisation de decoupe ss2 */
	for (i=0 ; i <= idx2 ; i++)
	{
		tab_decoupe[i]=1;
	}
	/* Elimination des zones de structure secondaire courte */
	for (i=0 ; i < j ; i=i+2)
	{
		//printf("%d %d\n",tab_u_ss2[i],tab_u_ss2[i+1]);
		for (k=tab_u_ss2[i] ; k < tab_u_ss2[i+1] ; k++)
		{
			tab_decoupe[k]=0;
		}
	}

	//Modification de la matrice de contact pour ne prendre en compte que les ss2

	if (ONLYSS2 == 1)
	{
		for (i=0 ; i <= idx2 ; i++)
		{
			for (j=0 ; j <= idx2 ;j++)
			{
				//printf("%d %d %d %d %f ",i,j,tab_ss2[i],tab_ss2[j], tab_pcontact[i][j] );
				if (tab_ss2[i] == 0 || tab_ss2[j] ==0)
				{
					tab_pcontact[i][j]=tab_pcontact[i][j]/10;
				}
				//printf("%f\n", tab_pcontact[i][j]);
			}
		}
	}

	//for (i=0 ; i <= idx2 ; i++)
	// {
	//	printf("%d %d\n",i, tab_decoupe[i]);
	// }
	//printf("dssp ok !\n");
}


/******************************************************************/
/*                                                                */
/*              Fonction Simple decoupage                         */
/*                                                                */
/******************************************************************/

void simple_cutting()
{
	int k,l;
	int i,j;
	int start_ju,i1,i2,end_ju; /* Valeur pour le rabottage */
	int ju1,ju2,ju3,ju4;
	float a,b,c;
	float coeff_matthews;
	int coo_max_i;
	int coo_max_j;
	int coo_min_i;
	int coo_min_j;
	int min_seg_size=MIN_SIZE_PU;
	int size_pu1,size_pu2;

	int i_p;

	/*Calcul de pcontact */
	pcontact=0.;
	for (k=start ; k <= end ; k++)
	{
		for (l=start ; l <= end ; l++)
		{
			pcontact=tab_pcontact[k][l]+pcontact;
		}
	}
	//printf("START:%d END:%d\n",start,end);
	/*Pour chaque position*/
	for (i=start ; i < end ; i++)
	{
		/*Si la position est decoupable on continu */
		i_p=i+1;
		if (tab_decoupe[i_p]==0) {continue;}
		/*Si taille des pu est superieure a min size on continue */
		size_pu1=i-start;
		size_pu2=end-i;
		if (size_pu1 >= min_seg_size && size_pu2 >= min_seg_size)
		{
			start_ju=start;
			i1=i;
			i2=i+1;
			end_ju=end;

			k=0;
			l=0;
			a=0.;
			b=0.;
			c=0.;

			/*Calcul du Coefficient de Matthews*/
			/* Calcul de "a"*/
			for (k=start_ju ; k <= i1 ; k++)
			{
				for (l=start_ju ; l <= i1 ; l++)
				{
					a=a+tab_pcontact[k][l];
				}
			}
			/* Calcul de "b"*/
			for (k=i2 ; k <= end_ju ; k++)
			{
				for (l=i2 ; l <= end ; l++)
				{
					b=b+tab_pcontact[k][l];
				}
			}

			/*Calcul de c*/
			/* c=(pcontact-a-b)/2; Si il n'y avait pas de ju */
			for (k=start_ju ; k <= i1 ; k++)
			{
				for (l=i2 ; l <=end_ju ; l++)
				{
					c=c+tab_pcontact[k][l];
				}
			}

			/*Calcul du coefficient de Mattews*/
			coeff_matthews=(a*b-c*c)/((a+c)*(b+c));
			//printf ("start:%d i1:%d i2:%d end:%d coeff:%f\n",start_ju,i1,i2,end_ju,coeff_matthews);
			/*Comptabilisation de la sone de coupure si le coeff est + eleve*/
			if (coeff_matthews > max_coeff_matthews)
			{
				max_start=start_ju;
				max_i1=i1;
				max_i2=i2;
				max_end=end_ju;
				max_coeff_matthews=coeff_matthews;
				nbre_de_coupe=1;
				best_pu=current_pu;
			}
		}
	}
	/*Imprimer coeff de matthews max et point de coupure idéal selon coeff*/
	//	printf("Simple decoupage: %d-%d %d-%d CM:%f\n",max_start,max_i1,max_i2,max_end,max_coeff_matthews);

}

/******************************************************************/
/*                                                                */
/*              Fonction Double decoupage                         */
/*                                                                */
/******************************************************************/

void double_cutting()
{

	int i,j;
	int k,l;
	int start_ju,i1,i2,j1,j2,end_ju; /* Valeur pour le rabottage */
	int ju1,ju2,ju3,ju4;
	float a,b,c;
	float coeff_matthews;
	int coo_max_i;
	int coo_max_j;
	int coo_min_i;
	int coo_min_j;
	int min_seg_size=MIN_SIZE_PU;
	int size_pu2;
	int i_p;
	int j_p;

	coo_max_i=end - min_seg_size;
	coo_max_j=end -(int)(min_seg_size/2);
	coo_min_i=start+min_seg_size-1;
	coo_min_j=0;

	for (i=coo_min_i ; i < coo_max_i ; i++)
	{
		/*Si la position est decoupable on continu */
		i_p=i+1;
		if (tab_decoupe[i_p]==0) {continue;}
		//printf("%d %d\n",i_p,tab_decoupe[i_p]);
		coo_min_j=i+min_seg_size;
		if (coo_min_j+min_seg_size < coo_max_j)
		{
			for (j=coo_min_j ; j <= coo_max_j ; j++)
			{
				/*Si la position est decoupable on continu */
				j_p=j+1;
				if (tab_decoupe[j_p]==0) {continue;}
				/*Si taille des pu est superieure a min size on continue */
				size_pu2=j-i;
				if(size_pu2 >= min_seg_size)
				{
					start_ju=start;
					i1=i;
					//if (tab_decoupe[i1]==0){continue;}
					i2=i+1;
					if (tab_decoupe[i1]==0){continue;}
					j1=j;
					//if (tab_decoupe[j1]==0){continue;}
					j2=j+1;
					if (tab_decoupe[j2]==0){continue;}
					end_ju=end;
					k=0;
					l=0;
					a=0;
					b=0;
					c=0;

					/*Calcul du Coefficient de Matthews*/
					/* Calcul de "a"*/
					for (k=i2 ; k <= j1 ; k++)
					{
						for (l=i2 ; l <= j1 ; l++)
						{
							a=a+tab_pcontact[k][l];
						}
					}
					/* Calcul de "b" 1*/
					for (k=start_ju ; k <= i1 ; k++)
					{
						for (l=start_ju ; l <= i1  ; l++)
						{
							b=b+tab_pcontact[k][l];
						}
					}
					/* Calcul de "b" 2*/
					for (k=j2 ; k <= end_ju ; k++)
					{
						for (l=j2 ; l <= end_ju ; l++)
						{
							b=b+tab_pcontact[k][l];
						}
					}

					/* Calcul de "b" 3 + "b" 4*/
					for (k=start_ju ; k <= i1 ; k++)
					{
						for (l=j2 ; l <= end_ju ; l++)
						{
							b=b+2*tab_pcontact[k][l];
						}
					}

					/*Calcul de c*/
					/*c=(pcontact-a-b)/2; */
					/*c1*/
					for (k=i2 ; k <= j1 ; k++)
					{
						for (l=start_ju ; l <= i1 ; l++)
						{
							c=c+tab_pcontact[k][l];
						}
					}
					/*c2*/
					for (k=j2 ; k <= end_ju ; k++)
					{
						for (l=i2 ; l <= j1 ; l++)
						{
							c=c+tab_pcontact[k][l];
						}
					}

					/*Calcul du coefficient de Mattews*/
					coeff_matthews=(a*b-c*c)/((a+c)*(b+c));
					//printf ("start:%d i1:%d i2:%d j1:%d j2:%d end:%d coeff:%f\n",start_ju,i1,i2,j1,j2,end_ju,coeff_matthews);
					/*Comptabilisation de la sone de coupure si le coeff est + eleve*/
					if (coeff_matthews > max_coeff_matthews)
					{
						max_start=start_ju;
						max_i1=i1;
						max_i2=i2;
						max_j1=j1;
						max_j2=j2;
						max_end=end_ju;
						max_coeff_matthews=coeff_matthews;
						nbre_de_coupe=2;
						best_pu=current_pu;
					}
				}
			}
		}
	}
	/*Imprimer coeff de matthews max et point de coupure ideal selon coeff*/
	//printf("double decoupage:%d-%d %d-%d %d-%d CM:%f best_pu:%d\n",max_start,max_i1,max_i2,max_j1,max_j2,max_end,max_coeff_matthews,best_pu);

}

/******************************************************************/
/*                                                                */
/*              Stockage des decoupes realisees                   */
/*                                                                */
/******************************************************************/

void save_pu()
{
	int i;
	int x;

	if (nbre_de_coupe==1)
	{
		/* Premiere PU */
		pu[iteration][new_nbre_pu][0]=max_start;
		pu[iteration][new_nbre_pu][1]=max_i1;
		/* Deuxieme PU */
		new_nbre_pu++;
		pu[iteration][new_nbre_pu][0]=max_i2;
		pu[iteration][new_nbre_pu][1]=max_end;
	}
	else if (nbre_de_coupe==2)
	{
		/* Premiere PU */
		pu[iteration][new_nbre_pu][0]=max_start;
		pu[iteration][new_nbre_pu][1]=max_i1;
		/* Deuxieme PU */
		new_nbre_pu++;
		pu[iteration][new_nbre_pu][0]=max_i2;
		pu[iteration][new_nbre_pu][1]=max_j1;
		/* Toisieme PU */
		new_nbre_pu++;
		pu[iteration][new_nbre_pu][0]=max_j2;
		pu[iteration][new_nbre_pu][1]=max_end;
	}

	/* Ajouter les old pu si il y lieu */
	if (iteration > 1)
	{
		for (x=0 ; x <= nbre_pu ; x++)
		{
			if (x == best_pu)
			{
			}
			else
			{
				new_nbre_pu++;
				pu[iteration][new_nbre_pu][0]=pu[iteration-1][x][0];
				pu[iteration][new_nbre_pu][1]=pu[iteration-1][x][1];
			}
		}

	}

	//for (x=0 ; x <= new_nbre_pu ; x++){printf("ITERATION:%d PU:%d %d-%d\n",iteration,x,pu[iteration][x][0],pu[iteration][x][1]);}




}

/******************************************************************/
/*                                                                */
/*              Fonction Calcul coeff homogeneitee                */
/*                                                                */
/******************************************************************/

double homogeneity(int start_pu , int end_pu)
{
	int k,l;
	int i,j;

	double pcontact1=0.0;
	double pcontact2=0.0;
	double pnormalize1=0.0;
	double pnormalize2=0.0;
	double H_1=0.0;
	double H_2=0.0;
	double Neq1=0.0;
	double Neq2=0.0;

	float seuil=0.5;

	int delta=0;

	int n=0;

	double I=0.0;

	/*Calcul de pcontact */

	/* Calcul probabilite de contact totale de la PU*/
	k=0;
	l=0;
	for (k=start_pu ; k <= end_pu ; k++)
	{
		for (l=start_pu ; l <= end_pu ; l++)
		{
			//if (k == l)continue;
			if (tab_pcontact[k][l] > seuil )
			{
				/*printf("proba:%f\n",tab_pcontact[k][l]);*/
				pcontact1=tab_pcontact[k][l]+pcontact1;
				delta=abs(k-l);
				if ( delta < 6)
				{
					pcontact2=tab_pcontact[k][l]+pcontact2;
				}

				/*printf("k:%d l:%d pcontact1 = %f tab_pcontact=%f\n",k,l,pcontact1, tab_pcontact[k][l]);*/
			}
		}
	}

	/* Normalisation (pcontact'=pcontact/pcontact_total */

	/*printf("START:%d END:%d\n",start,end);*/

	/*Pour chaque position*/
	H_1=0.;
	H_2=0.;
	for (i=start_pu ; i < end_pu ; i++)
	{
		for (j=start_pu ; j < end_pu ; j++)
		{
			//if (i == j)continue;
			if (tab_pcontact[i][j] <  0.0001) continue;
			if (tab_pcontact[i][j] > seuil )
			{
				/* Normalisation 1*/
				pnormalize1=(tab_pcontact[i][j]/pcontact1);
				pnormalize2=(tab_pcontact[i][j]/pcontact2);
				/*printf ("i:%d j:%d  PNORMALIZE1=%f PNORMALIZE2:%f\n",i,j,pnormalize1,pnormalize2);    */
				/* Calcul de l'entropie H_1 */
				H_1=H_1 + (pnormalize1 * log(pnormalize1));
				/*printf("LOG: %f\n",pnormalize1*log(pnormalize1));*/
				/* Calcul de l'entropie H_2 */
				delta=abs(i-j);
				if ( delta < 6)
				{
					H_2=H_2 + (pnormalize2 * log(pnormalize2));
				}
			}
		}
	}

	/* Calcul du coeff I d homogeneite */

	/* Calcul Neq (Nequivalent contact */
	//printf ("H_1: %f\n",H_1);
	//printf ("H_2: %f\n",H_2);
	Neq1=exp(-(H_1));
	Neq2=exp(-(H_2));
	//printf ("Neq1=%f Neq2=%f\n",Neq1,Neq2);
	/* Calcul I */
	n=end_pu-start_pu;
	I=(Neq1-Neq2)/n;
	//printf ("Homogeneity :%f  start:%d end:%d N n=%d\n",I,start_pu,end_pu,n);
	return(I);
}





/******************************************************************/
/*                                                                */
/*              Mutual information                                */
/*                                                                */
/******************************************************************/

void mutual_information()
{
	int x,y,x1,x2,y1,y2,k1,k2;
	float prob_zone[64][64];
	float sprob_zone[64];
	float sprob_tot;
	float entropie;
	float R;

	/* Balayage de toute la proteine et calcul des probabilites pour chaque pu vs chaque pu*/
	for (x=0 ; x <= new_nbre_pu ; x++)
	{
		x1=pu[iteration][x][0];
		x2=pu[iteration][x][1];
		for (y=0 ; y <= new_nbre_pu ; y++)
		{
			y1=pu[iteration][y][0];
			y2=pu[iteration][y][1];
			prob_zone[x][y]=0.;
			for (k1=x1 ; k1 <= x2 ; k1++)
			{
				for (k2=y1 ; k2 <= y2 ;k2++)
				{
					prob_zone[x][y]=prob_zone[x][y]+tab_pcontact[k1][k2];
				}
			}
		}

	}
	/* Calcul de Somme totale des probabilites et des sommes sur chaque ligne et colonne*/
	sprob_tot=0.;
	for (x=0 ; x <= new_nbre_pu ; x++)
	{
		sprob_zone[x]=0.;
		for (y=0 ; y <= new_nbre_pu ; y++)
		{
			sprob_zone[x]=sprob_zone[x]+prob_zone[x][y];
		}
		sprob_tot=sprob_tot+sprob_zone[x];
	}

	/* Calcul des rapports de probabilites*/
	for (x=0 ; x <= new_nbre_pu ; x++)
	{
		for (y=0 ; y <= new_nbre_pu ; y++)
		{
			prob_zone[x][y]=prob_zone[x][y] / sprob_tot;
		}
		sprob_zone[x]=sprob_zone[x]/sprob_tot;
	}

	entropie=0.;
	/* Calcul de l'information mutuelle */
	for (x=0 ; x <= new_nbre_pu ; x++)
	{
		for (y=0 ; y <= new_nbre_pu ; y++)
		{
			if (prob_zone[x][y] <= 0.00001 || sprob_zone[x] <= 0.00001 || sprob_zone[y] <= 0.00001)
			{
			}
			else
			{
				entropie=entropie+prob_zone[x][y]*log(prob_zone[x][y]/(sprob_zone[x]*sprob_zone[y]));
			}
		}
	}
	//CI=100*sqrt((1-exp(-2*entropie)));
	CI=100*sqrt(1-exp(-2*entropie));
	R=100*(1-exp(-2*entropie));
	printf("X.XX X.XX %f %f N ",CI,R);
	for (x=0 ; x <= new_nbre_pu ; x++){printf("%d %d ",pu[iteration][x][0]+1,pu[iteration][x][1]+1);}
	printf("\n");
	if (CI > MAXR2){exit (1);}

}
