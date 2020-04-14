//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/** @file cal_path.c
	@brief Calculate the path of where stuff is, and if it can't find a file look harder.  Win/Linux.
*/


#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "cal_path.h"
#include "util.h"
#include "inp.h"
#include <log.h>

#include <unistd.h>
#include <dirent.h>
#include <fcntl.h>
#include <stdarg.h>
#include <gpvdm_const.h>


#include <limits.h>
#include <sys/types.h>
#include <pwd.h>

int get_delta_path(struct simulation *sim,char *out, char *root,char *file_name)
{
int root_len=strlen(root);
int file_name_len=strlen(file_name);
if (root_len>file_name_len)
{
	strcpy(out,file_name);
	return -1;
}

if (root_len==0)
{
	strcpy(out,file_name);
	return 0;
}
if (strcmp_begin(file_name,root)==0)
{
	strcpy(out,file_name+root_len+1);
	return 0;
}

return -1;
}

int find_dll(struct simulation *sim, char *lib_path,char *lib_name)
{
char full_name[PATH_MAX];
char temp[PATH_MAX];
sprintf(full_name,"%s.so",lib_name);

join_path(2,lib_path,get_plugins_path(sim),full_name);
if (isfile(lib_path)==0)
{
	return 0;
}

struct dirent *next_file;
DIR *theFolder;

theFolder = opendir(get_plugins_path(sim));
if (theFolder!=NULL)
{
	while((next_file=readdir(theFolder))!=NULL)
	{
		split_dot(temp, next_file->d_name);
		if (strcmp(lib_name,temp)==0)
		{
			join_path(2,lib_path,get_plugins_path(sim),next_file->d_name);
			if (isfile(lib_path)==0)
			{
				closedir (theFolder);
				return 0;
			}

		}
	}

closedir (theFolder);

}

ewe(sim,"I can't find the dll %s,\n",lib_name);

return -1;
}

void set_path(struct simulation *sim,char *out, char *name)
{
char cwd[PATH_MAX];

char temp[PATH_MAX];



	//Check the home dir first
	join_path(3,temp,sim->home_path,"gpvdm_local",name);
	if ((isdir(temp)==0)||(isfile(temp)==0))
	{
		strcpy(out,temp);
		//printf(">>>>>>>>%s",temp);
		return;
	}

	//Check the cwd
	if (getcwd(cwd,PATH_MAX)==NULL)
	{
		ewe(sim,"cwd returned NULL, check if the directory exists.\n");
	}

	join_path(2,temp,cwd,name);

	if ((isdir(temp)==0)||(isfile(temp)==0))
	{
		strcpy(out,temp);
		return;
	}

	//check gpvdm_*
	join_path(3,temp,cwd,"gpvdm_data",name);

	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(3,temp,cwd,"gpvdm_core",name);

	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	//search the exe path
	join_path(2,temp,sim->exe_path,name);

	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	//search the exe path minus one level

	join_path(3,temp,sim->exe_path_dot_dot,"gpvdm_data",name);

	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(3,temp,sim->exe_path_dot_dot,"gpvdm_core",name);

	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(2,temp,"/usr/lib/gpvdm/",name);
	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(2,temp,"/usr/lib64/gpvdm/",name);
	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(2,temp,"/usr/share/gpvdm/",name);
	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	//Ubuntu
	join_path(2,temp,"/usr/lib/x86_64-linux-gnu/gpvdm/",name);
	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	join_path(2,temp,sim->share_path,name);
	if (isdir(temp)==0)
	{
		strcpy(out,temp);
		return;
	}

	if ((strcmp(name,"settings.inp")!=0)&&(strcmp(name,"info.inp")!=0))
	{
		ewe(sim,"I can't find the %s\n",name);
	}
}

void get_file_name_from_path(char *out,char *in)
{
	int i=0;

	strcpy(out,in);

	if (strlen(in)==0)
	{
		return;
	}

	for (i=strlen(in)-1;i>=0;i--)
	{
		if ((in[i]=='\\') || (in[i]=='/'))
		{
			i++;
			strcpy(out,(char*)(in+i));
			return;
		}
	}

}

int is_dir_in_path(char *long_path, char* search_dir)
{
	if( strstr(long_path, search_dir) != NULL)
	{
		return 0;
	}else
	{
		return -1;
	}
return -1;
}

void get_nth_dir_name_from_path(char *out,char *in,int n)
{
	int i=0;
	int ii=0;
	int pos=0;
	int start=0;
	int stop=0;
	int count=0;
	strcpy(out,in);

	if (strlen(in)==0)
	{
		return;
	}

	for (i=0;i<strlen(in);i++)
	{
		if ((in[i]=='\\') || (in[i]=='/') || (i==strlen(in)-1))
		{
			if (i!=0)
			{
				stop=i;
				if (i==strlen(in)-1)
				{
					stop++;
				}


				if (count==n)
				{

					for (ii=start;ii<stop;ii++)
					{
						out[pos]=in[ii];
						pos++;
					}
					out[pos]=0;
					return;
				}
				start=stop+1;
				count++;
			}else
			{
				start=i+1;		// move it past the first / in a unix string
			}
		}
	}

}
void cal_path(struct simulation *sim)
{
char cwd[PATH_MAX];
char temp[PATH_MAX];

strcpy(cwd,"");
strcpy(temp,"");

strcpy(sim->share_path,"nopath");

strcpy(sim->plugins_path,"");
strcpy(sim->lang_path,"");
strcpy(sim->input_path,"");
strcpy(sim->output_path,"");


memset(temp, 0, PATH_MAX * sizeof(char));
int len = readlink("/proc/self/exe", temp, PATH_MAX);
if (len == -1)
{
	ewe(sim,"IO error\n");
}


struct passwd *pw = getpwuid(getuid());

strcpy(sim->home_path,pw->pw_dir);


get_dir_name_from_path(sim->exe_path, temp);
get_dir_name_from_path(sim->exe_path_dot_dot, sim->exe_path);

if (isfile("configure.ac")==0)
{
	strcpy(sim->share_path,cwd);
	//printf_log(sim,"share path: %s\n",sim->share_path);
}else
if (isfile("ver.py")==0)
{
	path_up_level(temp, cwd);
	strcpy(sim->share_path,temp);
	//printf_log(sim,"share path: %s\n",sim->share_path);
}else
{
	strcpy(sim->share_path,"/usr/lib64/gpvdm/");
}

if (getcwd(cwd,PATH_MAX)==NULL)
{
	ewe(sim,"cwd returned NULL\n");
}

strcpy(sim->root_simulation_path,cwd);
strcpy(sim->output_path,cwd);
strcpy(sim->input_path,cwd);
set_path(sim,sim->plugins_path, "plugins");
//set_path(sim,sim->lang_path, "lang");
strcpy(sim->lang_path,"langdisabled");
set_path(sim,sim->materials_path, "materials");
set_path(sim,sim->cie_color_path, "cie_color");
set_path(sim,sim->shape_path, "shape");
set_path(sim,sim->emission_path, "emission");

set_path(sim,sim->spectra_path, "spectra");
//join_path(3,sim->cache_path,sim->home_path,"gpvdm_local","cache");
join_path(2,sim->cache_path,sim->input_path,"cache");
join_path(2,sim->gpvdm_local_path,sim->home_path,"gpvdm_local");

join_path(2,sim->tmp_path,sim->gpvdm_local_path,"tmp");
}



char *get_cache_path(struct simulation *sim)
{
return sim->cache_path;
}

char *get_gpvdm_local_path(struct simulation *sim)
{
return sim->gpvdm_local_path;
}

char *get_spectra_path(struct simulation *sim)
{
return sim->spectra_path;
}


char *get_materials_path(struct simulation *sim)
{
return sim->materials_path;
}

char *get_cie_color_path(struct simulation *sim)
{
return sim->cie_color_path;
}

char *get_shape_path(struct simulation *sim)
{
return sim->shape_path;
}

char *get_plugins_path(struct simulation *sim)
{
return sim->plugins_path;
}

char *get_lang_path(struct simulation *sim)
{
return sim->lang_path;
}

char *get_input_path(struct simulation *sim)
{
return sim->input_path;
}

char *get_output_path(struct simulation *sim)
{
return sim->output_path;
}

void set_output_path(struct simulation *sim,char *in)
{
strcpy(sim->output_path,in);
}

void set_input_path(struct simulation *sim,char *in)
{
strcpy(sim->input_path,in);
}

char *get_tmp_path(struct simulation *sim)
{
return sim->tmp_path;
}

void join_path(int max, ...)
{
	max=max+1;
	char temp[PATH_MAX];
	strcpy(temp,"");
	va_list arguments;
	int i;
	va_start ( arguments, max );
	char *ret=va_arg ( arguments, char * );
	strcpy(ret,"");
	for (i = 1; i < max; i++ )
	{
		if ((i!=1)&&(strcmp(temp,"")!=0))
		{
			strcat(ret,"/");
		}
		strcpy(temp,va_arg ( arguments, char * ));
		strcat(ret,temp);
	}
	va_end ( arguments );                  // Cleans up the list

	return;
}


/**Make sure the slashes go the right way in a string for which ever OS we are on.
@param path path to check
*/
void assert_platform_path(char * path)
{
	int i=0;
	char temp[PATH_MAX];
	strcpy(temp,"");
	int max=strlen(path);
	for (i=0;i<max;i++)
	{
		if ((path[i]=='\\')||(path[i]=='/'))
		{
			strcat(temp,"/");
		}else
		{
			temp[i]=path[i];
			temp[i+1]=0;
		}


	}

	strcpy(path,temp);

	return;
}
