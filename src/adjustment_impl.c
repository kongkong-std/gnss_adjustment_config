#include "../include/adjustment_impl.h"

void ImplNetworkAdjustment(Config *var_conf, AdjGraph *var_data)
{
    if (var_conf->mtd_Weight == 0)
    {
        // equal weight
        puts(">>>>>>>> equal weight adjustment:");
        ImplNetworkAdjustment_0(var_conf, var_data);
    }
    else if (var_conf->mtd_Weight == 1)
    {
        // diagonal weight
        puts(">>>>>>>> diagonal weight adjustment:");
        ImplNetworkAdjustment_1(var_conf, var_data);
    }
    else if (var_conf->mtd_Weight == 2)
    {
        // full weight
        puts(">>>>>>>> full weight adjustment:");
        ImplNetworkAdjustment_2(var_conf, var_data);
    }
}

void ImplNetworkAdjustment_0(Config *var_conf, AdjGraph *var_data)
{
}

void ImplNetworkAdjustment_1(Config *var_conf, AdjGraph *var_data)
{
}

void ImplNetworkAdjustment_2(Config *var_conf, AdjGraph *var_data)
{
}