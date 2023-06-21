/**
 * @file main_function.c
 * @author Zikang Qin
 * @brief function in main function,
 * mainly involves configure file process and
 * source data process
 *
 * @version 0.1
 * @date 2023-06-21
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "../include/main.h"

/**
 * @callgraph
 * @brief source data process, assigning values
 * to graph data structure with source data file
 *
 * @param [in] fp_data FILE pointer
 * @param [in] file_data path to source data file
 * @param [in] var_conf configure data for gnss adjustment of network
 * @param [in,out] var_data graph data structure
 */
void ProcessSourceData(FILE *fp_data, char *file_data, Config *var_conf, AdjGraph *var_data)
{
	if ((fp_data = fopen(file_data, "rb")) == NULL)
	{
		fprintf(stderr, "Cannot open file \"%s\"\n", file_data);
		exit(EXIT_FAILURE);
	}

	int node_location[2] = {0};	 // node enumeration
	double weight_data[9] = {0}; // weight data in edge
	for (int index = 0; index < var_conf->cnt_BaseLine; index++)
	{
		for (int index_i = 0; index_i < 2; index_i++)
		{
			fscanf(fp_data, "%d", node_location + index_i);
		}

		EncodeStationId(var_conf, 2, node_location);

		for (int index_i = 0; index_i < 9; index_i++)
		{
			fscanf(fp_data, "%lg", weight_data + index_i);
		}
		(var_data->count_edge)++; // count of edge augment by 1
		GraphInsert(var_data, node_location, weight_data, 9);
	}

	fclose(fp_data);

#if 1 // graph display
	puts(">>>>>>>>>>>>graph display:");
	GraphDisplay(var_data);
#endif
}

/**
 * @callgraph
 * @brief encoding base station and rover station
 * sequentially start from 0
 *
 * @param [in] var_conf configure data for gnss adjustment of network
 * @param [in] len_node count of stations
 * @param [in,out] node array of station code
 */
void EncodeStationId(Config *var_conf, int len_node, int *node)
{
	for (int index = 0; index < len_node; index++)
	{
		for (int index_i = 0; index_i < var_conf->cnt_BaseStation; index_i++)
		{
			if (*(node + index) == var_conf->id_BaseStation[index_i])
			{
				*(node + index) = index_i;
			}
		}
		for (int index_i = 0; index_i < var_conf->cnt_RoverStation; index_i++)
		{
			if (*(node + index) == var_conf->id_RoverStation[index_i])
			{
				*(node + index) = index_i + var_conf->cnt_BaseStation;
			}
		}
	}
}

/**
 * @callgraph
 * @brief free memory of configure struct
 *
 * @param [in,out] adjust_conf configure data for gnss adjustment of network
 */
void free_conf(Config *adjust_conf)
{
#if 0
	if (adjust_conf->id_BaseStation != NULL)
	{
		free(adjust_conf->id_BaseStation);
	}
	if (adjust_conf->coo_BaseStation != NULL)
	{
		free(adjust_conf->coo_BaseStation);
	}
	if (adjust_conf->id_RoverStation != NULL)
	{
		free(adjust_conf->id_RoverStation);
	}
#endif
	free(adjust_conf);
}

/**
 * @callgraph
 * @brief process sorted configure file,
 * this type of file contains only data
 * information, such as, count of stations,
 * base station coordinate and etc.
 * assigning values to configure struct with
 * configure file
 *
 * @param [in] fp_conf FILE pointer
 * @param [in] file_conf path of configure file
 * @param [in,out] var_conf configure struct
 */
void ProcessConfigureFile_sort(FILE *fp_conf, char *file_conf, Config *var_conf)
{
	if ((fp_conf = fopen(file_conf, "rb")) == NULL)
	{
		fprintf(stderr, "Cannot open file \"%s\"\n", file_conf);
		exit(EXIT_FAILURE);
	}

	fscanf(fp_conf, "%d", &var_conf->cnt_Station);
	fscanf(fp_conf, "%d", &var_conf->cnt_BaseStation);
	fscanf(fp_conf, "%d", &var_conf->cnt_RoverStation);

	for (int index = 0; index < var_conf->cnt_BaseStation; index++)
	{
		fscanf(fp_conf, "%d", var_conf->id_BaseStation + index);
	}
	for (int index = 0; index < var_conf->cnt_RoverStation; index++)
	{
		fscanf(fp_conf, "%d", var_conf->id_RoverStation + index);
	}
	for (int index = 0; index < 3 * var_conf->cnt_BaseStation; index++)
	{
		fscanf(fp_conf, "%lg", var_conf->coo_BaseStation + index);
	}

	fscanf(fp_conf, "%d", &var_conf->cnt_BaseLine);
	fscanf(fp_conf, "%d", &var_conf->mtd_Weight);

	fclose(fp_conf);
}

/**
 * @callgraph
 * @brief process configure file,
 * assigning values to configure struct
 * with configure file
 *
 * @param [in] fp_conf FILE pointer
 * @param [in] file_conf path of configure file
 * @param [in,out] var_conf configure struct
 */
void ProcessConfigureFile(FILE *fp_conf, char *file_conf, Config *var_conf)
{
	if ((fp_conf = fopen(file_conf, "rb")) == NULL)
	{
		fprintf(stderr, "Cannot open file \"%s\"!\n", file_conf);
		exit(EXIT_FAILURE);
	}

	int cnt_id_BaseStation = 0, cnt_id_RoverStation = 0, cnt_coo_BaseStation = 0;

	while (!feof(fp_conf))
	{
		char buffer[MAXBUFFER];
		fgets(buffer, MAXBUFFER, fp_conf);
		if (strstr(buffer, "cnt_Station"))
		{
			// get cnt_Station
			sscanf(buffer, "%d", &var_conf->cnt_Station);
		}
		if (strstr(buffer, "cnt_BaseStation"))
		{
			// get cnt_BaseStation
			sscanf(buffer, "%d", &var_conf->cnt_BaseStation);
		}
		if (strstr(buffer, "cnt_RoverStation"))
		{
			// get cnt_RoverStation
			sscanf(buffer, "%d", &var_conf->cnt_RoverStation);
		}
		if (strstr(buffer, "id_BaseStation"))
		{
			// get id_BaseStation
			sscanf(buffer, "%d", var_conf->id_BaseStation + cnt_id_BaseStation);
			cnt_id_BaseStation++;
		}
		if (strstr(buffer, "id_RoverStation"))
		{
			// get id_RoverStation
			sscanf(buffer, "%d", var_conf->id_RoverStation + cnt_id_RoverStation);
			cnt_id_RoverStation++;
		}
		if (strstr(buffer, "coo_BaseStation"))
		{
			// get coo_BaseStation
			sscanf(buffer, "%lf", var_conf->coo_BaseStation + cnt_coo_BaseStation);
			cnt_coo_BaseStation++;
		}
		if (strstr(buffer, "cnt_BaseLine"))
		{
			// get cnt_BaseLine
			sscanf(buffer, "%d", &var_conf->cnt_BaseLine);
		}
		if (strstr(buffer, "mtd_Weight"))
		{
			// get mtd_Weight
			sscanf(buffer, "%d", &var_conf->mtd_Weight);
		}
	}

	fclose(fp_conf);
}
