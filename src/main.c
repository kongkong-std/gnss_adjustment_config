#include "../include/main.h"

int main(int argc, char **argv)
{
    char *file_configure = argv[1]; // configure file
    char *file_data = argv[2];      // data file

    Config *adjust_conf = NULL;
    adjust_conf = (Config *)malloc(sizeof(Config));
    memset(adjust_conf, 0, sizeof(Config));

    FILE *fp_adjust = NULL;

    // processing configure file
    ProcessConfigureFile(fp_adjust, file_configure, adjust_conf);
    // ProcessConfigureFile_sort( fp_adjust, file_configure, adjust_conf );

    puts("============ GNSS network generation ============");

#if 0
    printf("adjust_conf->cnt_Station = %d\n", adjust_conf->cnt_Station);
    printf("adjust_conf->cnt_BaseStation = %d\n", adjust_conf->cnt_BaseStation);
    printf("adjust_conf->cnt_RoverStation = %d\n", adjust_conf->cnt_RoverStation);
    for (int index = 0; index < adjust_conf->cnt_BaseStation; index++)
    {
        printf("adjust_conf->id_BaseStation[ %d ] = %d\n", index, *(adjust_conf->id_BaseStation + index));
    }
    for (int index = 0; index < adjust_conf->cnt_RoverStation; index++)
    {
        printf("adjust_conf->id_RoverStation[ %d ] = %d\n", index, *(adjust_conf->id_RoverStation + index));
    }
    for (int index = 0; index < adjust_conf->cnt_BaseStation * 3; index++)
    {
        printf("%.4lf ", adjust_conf->coo_BaseStation[index]);
        if ((index + 1) % 3 == 0)
        {
            putchar('\n');
        }
    }
    printf("adjust_conf->cnt_BaseLint = %d\n", adjust_conf->cnt_BaseLine);
    printf("adjust_conf->mtd_Weight = %d\n", adjust_conf->mtd_Weight);
#endif

    // reading source file
    AdjGraph *value_graph = GraphGeneration(adjust_conf->cnt_Station); // graph generation with count vertices
    // AdjGraph * value_graph = GraphGeneration( 6 );    // graph generation with count vertices
    ProcessSourceData(fp_adjust, file_data, adjust_conf, value_graph);

    // gnss network adjustment implmentation
    puts( "\n============ Implementation of GNSS network adjustment ============" );
    ImplNetworkAdjustment( adjust_conf, value_graph );

    // free memory
    free_conf(adjust_conf);
    GraphDestroy(value_graph);

    return 0;
}
