//
// General-purpose Photovoltaic Device Model gpvdm.com- a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
//
// Copyright (C) 2012-2017 Roderick C. I. MacKenzie info at gpvdm dot com
//
// https://www.gpvdm.com
//
//
// This program is free software; you can redistribute it and/or modify it
// under the terms and conditions of the GNU Lesser General Public License,
// version 2.1, as published by the Free Software Foundation.
//
// This program is distributed in the hope it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
//

/** @file util.c
	@brief Utility functions.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <fcntl.h>
#include "util.h"
#include "log.h"
#include <gpvdm_const.h>
#include <lang.h>
#include <math.h>
#include <ctype.h>
#include <cal_path.h>

static char* unused_pchar __attribute__((unused));

void gpvdm_mkdir(char *file_name)
{
struct stat st = {0};

	if (stat(file_name, &st) == -1)
	{
			mkdir(file_name, 0700);
	}

}




/**Get length of a file in lines
@param file_name file name
*/
int get_file_len(struct simulation *sim,char *file_name)
{
FILE *file;
if (!(file=fopen(file_name,"r")))
{
	printf_log(sim,"Error opening file %s\n",file_name);
	exit(0);
}
char buffer[1000];

int i;
i=0;
char *p;
do
{
	buffer[0]=0;
	unused_pchar=fgets(buffer, 1000, file);
	p=buffer;
	if (buffer[0]=='#')
	{
		i--;
	}else
	{
		//Check for empty line
		while(*p==' ' || *p=='\t') p++;
		if ((*p=='\r')||(*p=='\n')||(*p==0)) i--;
	}
i++;
}while(!feof(file));
//i--;
fclose(file);
return i;
}






void print_hex(struct simulation *sim,unsigned char *data)
{
int i;
for (i=0;i<16;i++)
{
	printf_log(sim,"%02x",data[i]);
}
printf_log(sim,"\n");
}


void write_x_y_to_file(struct simulation *sim,char *name,double *x,double *y,int len)
{
int i;
FILE *out;

	out=fopen(name,"w");
	if (out==NULL)
	{
		ewe(sim,"%s %s\n",_("Error writing file"),name);
	}

	for (i=0;i<len;i++)
	{
		fprintf(out,"%le %le\n",x[i],y[i]);
	}
	fclose(out);
}

void write_x_y_z_to_file(struct simulation *sim,char *name,double *x,double *y,double *z,int len)
{
int i;
FILE *out;

	out=fopen(name,"w");
	if (out==NULL)
	{
		ewe(sim,"%s %s\n",_("Error writing file"),name);
	}

	for (i=0;i<len;i++)
	{
		fprintf(out,"%le %le %le\n",x[i],y[i],z[i]);
	}
	fclose(out);
}

void str_to_lower(char *out, char *in)
{
	int i=0;
	for (i=0;i<strlen(in);i++)
	{
		out[i]=tolower(in[i]);
	}
	out[i]=0;

}

int check_int(char *in)
{
int i=0;
int numeric=TRUE;
for (i=0;i<strlen(in);i++)
{
	if ((in[i]<48)||(in[i]>57))
	{
		numeric=FALSE;
		break;
	}
}
return numeric;
}
static int unused __attribute__((unused));

int english_to_bin(struct simulation *sim, char * in)
{
char temp[100];
int ret=0;

str_to_lower(temp, in);

if (check_int(temp)==TRUE)
{
sscanf(temp,"%d",&ret);
return ret;
}

if (strcmp(temp,"true")==0)
{
	return TRUE;
}else
if (strcmp(temp,"false")==0)
{
	return FALSE;
}else
if (strcmp(temp,"1")==0)
{
	return TRUE;
}else
if (strcmp(temp,"0")==0)
{
	return FALSE;
}else
if (strcmp(temp,"yes")==0)
{
	return TRUE;
}else
if (strcmp(temp,"no")==0)
{
	return FALSE;
}else
if (strcmp(temp,"left")==0)
{
	return LEFT;
}else
if (strcmp(temp,"links")==0)
{
	return LEFT;
}else
if (strcmp(temp,"ja")==0)
{
	return TRUE;
}else
if (strcmp(temp,"nein")==0)
{
	return FALSE;
}else
if (strcmp(temp,"right")==0)
{
	return RIGHT;
}else
if (strcmp(temp,"rechts")==0)
{
	return RIGHT;

}else
if (strcmp(temp,"gaus")==0)
{
	return 0;
}else
if (strcmp(temp,"exp")==0)
{
	return 1;
}else
if (strcmp(temp,"exponential")==0)
{
	return dos_exp;
}else
if (strcmp(temp,"complex")==0)
{
	return dos_an;
}
else
if (strcmp(temp,"open_circuit")==0)
{
	return OPEN_CIRCUIT;
}else
if (strcmp(temp,"load")==0)
{
	return LOAD;
}else
if (strcmp(temp,"ideal_diode_ideal_load")==0)
{
	return IDEAL_DIODE_IDEAL_LOAD;
}else
if (strcmp(temp,"none")==0)
{
	return log_level_none;
}else
if (strcmp(temp,"screen")==0)
{
	return log_level_screen;
}else
if (strcmp(temp,"disk")==0)
{
	return log_level_disk;
}else
if (strcmp(temp,"screen_and_disk")==0)
{
	return log_level_screen_and_disk;
}else
if (strcmp(temp,"newton")==0)
{
	return FIT_NEWTON;
}else
if (strcmp(temp,"simplex")==0)
{
	return FIT_SIMPLEX;
}else
if (strcmp(temp,"bfgs")==0)
{
	return FIT_BFGS;
}else
if (strcmp(temp,"top")==0)
{
	return TOP;
}else
if (strcmp(temp,"bottom")==0)
{
	return BOTTOM;
}else
if (strcmp(temp,"right")==0)
{
	return RIGHT;
}else
if (strcmp(temp,"left")==0)
{
	return LEFT;
}else

if (strcmp(temp,"electron")==0)
{
	return ELECTRON;
}else
if (strcmp(temp,"hole")==0)
{
	return HOLE;
}else
if (strcmp(temp,"mb_equation")==0)
{
	return mb_equation;
}else
if (strcmp(temp,"mb_look_up_table")==0)
{
	return mb_look_up_table;
}else
if (strcmp(temp,"fd_look_up_table")==0)
{
	return fd_look_up_table;
}else
if (strcmp(temp,"mb_look_up_table_analytic")==0)
{
	return mb_look_up_table_analytic;
}else
if (strcmp(temp,"ohmic")==0)
{
	return contact_ohmic;
}else
if (strcmp(temp,"schottky")==0)
{
	return contact_schottky;
}else
if (strcmp(temp,"ray_run_never")==0)
{
	return ray_run_never;
}else
if (strcmp(temp,"ray_run_once")==0)
{
	return ray_run_once;
}else
if (strcmp(temp,"ray_run_step")==0)
{
	return ray_run_step;
}else
if (strcmp(temp,"ray_emission_single_point")==0)
{
	return ray_emission_single_point;
}else
if (strcmp(temp,"ray_emission_electrical_mesh")==0)
{
	return ray_emission_electrical_mesh;
}else
if (strcmp(temp,"dump_verbosity_everything")==0)
{
	return dump_verbosity_everything;
}else
if (strcmp(temp,"dump_verbosity_key_results")==0)
{
	return dump_verbosity_key_results;
}else
if (strcmp(temp,"measure_voltage")==0)
{
	return measure_voltage;
}else
if (strcmp(temp,"measure_current")==0)
{
	return measure_current;
}else
if (strcmp(temp,"contact")==0)
{
	return LAYER_CONTACT;
}else
if (strcmp(temp,"active")==0)
{
	return LAYER_ACTIVE;
}else
if (strcmp(temp,"other")==0)
{
	return LAYER_OTHER;
}else
if (strcmp(temp,"thermal_hydrodynamic")==0)
{
	return THERMAL_HYDRODYNAMIC;
}else
if (strcmp(temp,"thermal_simple")==0)
{
	return THERMAL_SIMPLE;
}


ewe(sim,"%s %s\n",_("I don't understand the command"),in);
return 0;
}

double read_value(struct simulation *sim,char *file,int skip,int line)
{
FILE *in;
char buf[1000];
double value;
in=fopen(file,"r");
if (in==NULL)
{
ewe(sim,"%s %s\n",_("Can not read file"),file);
}
int l=0;

do
{
l++;
	unused=fscanf(in,"%s",buf);

	if (l==line)
	{
		sscanf((buf+skip),"%le\n",&value);
		break;
	}


}while (!feof(in));

fclose(in);


return value;
}

void safe_file(struct simulation *sim,char *name)
{
FILE *file;
file = fopen(name, "rb");

if (!file)
{
	ewe(sim,"%s: %s\n",_("File not found"),name);
}

fclose(file);
}




void randomprint(struct simulation *sim,char *in)
{
	int i;

	wchar_t wide[1000];
	//char temp[1000];
	int len=mbstowcs(wide, in, 1000);
	//wprintf(L"%S",wide);
	for (i=0;i<len;i++)
	{
	int rnd=(float)5.0*rand()/(float)RAND_MAX;
		if (rnd==0) textcolor(sim,fg_wight);
		if (rnd==1) textcolor(sim,fg_red);
		if (rnd==2) textcolor(sim,fg_green);
		if (rnd==3) textcolor(sim,fg_yellow);
		if (rnd==4) textcolor(sim,fg_blue);
		if (rnd==5) textcolor(sim,fg_purple);

		if ((wide[i]!='\n')||(sim->html==FALSE))
		{
			//mbstowcs(wide, in, 1000);
			//swprintf(temp,L"%C",wide[i]);
			printf_log(sim,"%C",wide[i]);
		}else
		{
			printf_log(sim,"<br>");
		}

		textcolor(sim,fg_reset);

		}

fflush(stdout);
}

FILE *fopena(char *path,char *name,const char *mode)
{
char wholename[200];
join_path(2, wholename,path,name);

FILE *pointer;
pointer=fopen(wholename,mode);

return pointer;
}


int scanarg( char* in[],int count,char * find)
{
int i;
for (i=0;i<count;i++)
{
if (strcmp(in[i],find)==0) return TRUE;
}
return FALSE;
}

int get_arg_plusone_pos( char* in[],int count,char * find)
{
int i;
for (i=0;i<count;i++)
{
if (strcmp(in[i],find)==0)
{
       if ((i+1)<count)
       {
               return i+1;
       }else
       {
               return FALSE;
       }
}
}
return FALSE;
}

char * get_arg_plusone( char* in[],int count,char * find)
{
int i;
static char no[] = "";
for (i=0;i<count;i++)
{

if (strcmp(in[i],find)==0)
{
       if ((i+1)<count)
       {
               return in[i+1];
       }else
       {
               return no;
       }
}
}

return no;
}



void copy_file(struct simulation *sim,char *output,char *input)
{
char buf[8192];
struct stat results;
int in_fd = open(input, O_RDONLY);
if (in_fd== -1)
{
	ewe(sim,"%s: %s\n",_("Can not open file"),input);
}

stat(input, &results);

int out_fd = open(output, O_WRONLY | O_CREAT| O_TRUNC,results.st_mode);
if (in_fd== -1)
{
	ewe(sim,"%s: %s\n",_("Can not open file"),output);
}


while (1)
{
    ssize_t result = read(in_fd, buf, 8192*sizeof(char));

    if (result==0)
	{
		break;
	}
    write(out_fd, buf, result*sizeof(char));
}

close(in_fd);
close(out_fd);
}



int path_up_level(char *out, char *in)
{
int i=0;
strcpy(out,in);
int len=strlen(out);
if (len<1)
{
	return -1;
}

if (len!=3)
{
	if (out[len-1]=='\\')
	{
		out[len-1]=0;
		len=strlen(out);
	}
}

if (len!=1)
{
	if (out[len-1]=='/')
	{
		out[len-1]=0;
		len=strlen(out);
	}
}

for (i=len;i>=0;i--)
{

		if (out[i]=='\\')
		{
			out[i+1]=0;
			break;
		}



		if (out[i]=='/')
		{
			out[i+1]=0;
			break;
		}
}

return 0;
}


int get_dir_name_from_path(char *out, char *in)
{
strcpy(out,in);

int i=0;
int len=strlen(in);
for (i=len;i>0;i--)
{
	//#ifdef windows
		if ((in[i]=='\\')||(in[i]=='/'))
		{
			out[i]=0;
			return 0;
		}
	//#else
	//	if (in[i]=='/')
	//	{
	//		out[i]=0;
	//		return 0;
	//	}
	//#endif
}

if (len>0)
{
	out[0]=0;
	return 0;
}

return -1;
}

/**
 * @brief Tests for a directory, returns 0 if found non zero if not.
 *
 * Detailed explanation.
 */
int isdir(const char *path)
{
	struct stat statbuf;
	if (stat(path, &statbuf) != 0)
	{
		return 1;
	}else
	{
		if (S_ISDIR(statbuf.st_mode)==0)
		{
			return 1;
		}else
		{
			return 0;
		}
	}
}

/**
 * @brief Tests for a file, returns 0 if found non zero if not.
 *
 * Detailed explanation.
 */
int isfile(char *in)
{
	struct stat statbuf;
	if (stat(in, &statbuf) != 0)
	{
		return 1;
	}else
	{
		if (S_ISREG(statbuf.st_mode)==0)
		{
			return 1;
		}else
		{
			return 0;
		}
	}

}

void remove_file(struct simulation *sim,char* file_name)
{
	if (isfile(file_name)==0)
	{
		remove(file_name);
	}
}

void remove_dir_ittr(struct simulation *sim,char* dir_name,int depth)
{
	depth++;
	struct dirent *next_file;
	DIR *theFolder;
	char filepath[PATH_MAX];

	theFolder = opendir(dir_name);
	if (theFolder!=NULL)
	{
		while((next_file=readdir(theFolder))!=NULL)
		{
			if ((strcmp(next_file->d_name,".")!=0)&&(strcmp(next_file->d_name,"..")!=0))
			{
				join_path(2, filepath,dir_name,next_file->d_name);
				if (isdir(filepath)==0)
				{
					remove_dir_ittr(sim,filepath,depth);
					printf_log(sim,"%s =%s\n",_("Deleting directory"),filepath);
						remove(filepath);
				}else
				{
					remove_file(sim,filepath);
				}
			}
		}

		closedir (theFolder);

		if (depth==0)
		{
				remove(dir_name);
		}
	}
}

void remove_dir(struct simulation *sim,char* dir_name)
{
	remove_dir_ittr(sim,dir_name,-1);

}

/**This is a version of the standard fgets, but it will also accept a 0x0d as a new line.
@param buf output buffer
@param len max length of buffer
@param file file handle
*/
int gpvdm_fgets(char *buf,int len,FILE *file)
{
	char dat;
	int pos=0;

	if (feof(file))
	{
		return -1;
	}

	do
	{
		dat=(char)fgetc(file);
		if (feof(file))
		{
			break;
		}

		if ((dat=='\n')||(dat=='\r')||(dat==0x0d))
		{
			break;
		}

		buf[pos]=dat;

		pos++;

		if (pos>len)
		{
			break;
		}

	}while(1);
	buf[pos]=0;

	return pos;
}
