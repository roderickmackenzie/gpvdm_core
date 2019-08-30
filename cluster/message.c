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

/** @file message.c
@brief message passing functions
*/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <arpa/inet.h>
#include "inp.h"
#include "util.h"
#include <sys/ioctl.h>
#include <net/if.h>
#include "inp.h"
#include "tx_packet.h"

int send_message(char *message)
{
	char buf[LENGTH];

	struct node_struct* master=NULL;
	master=node_find_master();

	if (master!=NULL)
	{

		struct tx_struct packet;
		tx_struct_init(&packet);
		tx_set_id(&packet,"gpvdm_message");
		printf("%s\n",message);
		strncpy(packet.message,message,TX_STRUCT_DATA_ITEM_SIZE);
		packet.message[TX_STRUCT_DATA_ITEM_SIZE-1]=0;
		tx_packet(master->sock,&packet,buf);

	}else
	{
		return -1;
	}

return 0;
}



