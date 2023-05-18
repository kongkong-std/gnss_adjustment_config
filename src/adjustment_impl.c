#include "../include/adjustment_impl.h"

void LSEUpdateSolution(Config *var_conf, double *lse_sol, double *lse_b,
                       int lse_size_row, int lse_size_column, Solution *var_solution)
{
    for (int index = 0; index < var_conf->cnt_RoverStation; index++)
    {
        var_solution->id_RoverStation[index] = var_conf->id_RoverStation[index];
    }
    for (int index = 0; index < var_conf->cnt_RoverStation * 3; index++)
    {
        var_solution->coo_RoverStation[index] = lse_sol[index];
    }
    var_solution->norm_result = lse_norm_2( lse_b + lse_size_column, lse_size_row - lse_size_column );
}

double lse_norm_2(double *a, int dimension)
{
    double value = 0;

    for (int index = 0; index < dimension; index++)
    {
        value += a[index] * a[index];
    }

    return sqrt(value);
}

void ImplNetworkAdjustment(Config *var_conf, AdjGraph *var_data, Solution *var_solution)
{
    if (var_conf->mtd_Weight == 0)
    {
        // equal weight
        puts(">>>>>>>> equal weight adjustment:");
        ImplNetworkAdjustment_0(var_conf, var_data, var_solution);
    }
    else if (var_conf->mtd_Weight == 1)
    {
        // diagonal weight
        puts(">>>>>>>> diagonal weight adjustment:");
        ImplNetworkAdjustment_1(var_conf, var_data, var_solution);
    }
    else if (var_conf->mtd_Weight == 2)
    {
        // full weight
        puts(">>>>>>>> full weight adjustment:");
        ImplNetworkAdjustment_2(var_conf, var_data, var_solution);
    }
}

void ImplNetworkAdjustment_0(Config *var_conf, AdjGraph *var_data, Solution *var_solution)
{
    // data size
    int data_row = 0, data_column = 0;
    data_row = var_data->count_edge;
    data_column = sizeof(var_data->array->head->weight) / sizeof(var_data->array->head->weight[0]);

    int **vertex_enum = NULL;      // vertex enumeration of baseline data
    double **data_BaseLine = NULL; // baselien data, contains baseline vector, standard deviation

    /*
     * size vertex_enum = data_row x 2
     * size data_BaseLine = data_row x data_column
     */
    vertex_enum = (int **)malloc(data_row * sizeof(int *));
    for (int index = 0; index < data_row; index++)
    {
        *(vertex_enum + index) = (int *)malloc(2 * sizeof(int));
    }
    data_BaseLine = (double **)malloc(data_row * sizeof(double *));
    for (int index = 0; index < data_row; index++)
    {
        *(data_BaseLine + index) = (double *)malloc(data_column * sizeof(double));
    }

    // graph traversal
    for (int index = 0, index_temp = 0; index < var_data->count_vertex && index_temp < data_row; index++)
    {
        AdjListNode *head_node = var_data->array[index].head;

        while (head_node != NULL)
        {
            // vertex enumeration
            vertex_enum[index_temp][0] = index;
            vertex_enum[index_temp][1] = head_node->dest;

            // baseline data
            for (int index_i = 0; index_i < data_column; index_i++)
            {
                data_BaseLine[index_temp][index_i] = head_node->weight[index_i];
            }

            head_node = head_node->next;
            index_temp++;
        }
    }

#if 0 // print data
    puts(">>>>>>>>>>>>vertex_enum:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 2; index_j++)
        {
            printf("%d ", vertex_enum[index_i][index_j]);
        }
        putchar('\n');
    }
    puts("\n>>>>>>>>>>>>data_BaseLine:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            printf("%.4lf ", data_BaseLine[index_i][index_j]);
        }
        for( int index_j = 3; index_j < data_column; index_j++ )
        {
            printf( "%.4le ", data_BaseLine[ index_i ][ index_j ] );
        }
        putchar('\n');
    }
#endif

    // least square equation
    /*
     * sol = \arg\min || Ax - b ||
     * size A = data_row * 3 x cnt_RoverStation * 3
     * size sol = size x = cnt_RoverStation * 3 x 1
     * size b = data_row * 3 x 1
     * residual = Ax - b
     */
    double *lse_sol = NULL, *lse_b = NULL, *lse_residual = NULL;
    double **lse_A = NULL;
    int lse_size_column = var_conf->cnt_RoverStation * 3;
    int lse_size_row = data_row * 3;

    lse_sol = (double *)malloc(lse_size_column * sizeof(double));
    lse_b = (double *)malloc(lse_size_row * sizeof(double));
    lse_residual = (double *)malloc(lse_size_row * sizeof(double));
    lse_A = (double **)malloc(lse_size_row * sizeof(double *));
    for (int index = 0; index < lse_size_row; index++)
    {
        *(lse_A + index) = (double *)malloc(lse_size_column * sizeof(double));
    }

    InitLSE(lse_sol, lse_b, lse_residual, lse_A, lse_size_row, lse_size_column);

    LSESolver_0(var_conf, vertex_enum, data_BaseLine, data_row, data_column, 2,
                lse_A, lse_b, lse_size_row, lse_size_column, lse_sol, lse_residual);

    // updating solution
    LSEUpdateSolution(var_conf, lse_sol, lse_b, lse_size_row, lse_size_column, var_solution);

    // free memory
    for (int index = 0; index < lse_size_row; index++)
    {
        free(*(lse_A + index));
    }
    free(lse_A);
    free(lse_sol);
    free(lse_b);
    free(lse_residual);
    for (int index = 0; index < data_row; index++)
    {
        free(*(data_BaseLine + index));
    }
    free(data_BaseLine);
    for (int index = 0; index < data_row; index++)
    {
        free(*(vertex_enum + index));
    }
    free(vertex_enum);
}

void ImplNetworkAdjustment_1(Config *var_conf, AdjGraph *var_data, Solution *var_solution)
{
    // data size
    int data_row = 0, data_column = 0;
    data_row = var_data->count_edge;
    data_column = sizeof(var_data->array->head->weight) / sizeof(var_data->array->head->weight[0]);

    int **vertex_enum = NULL;      // vertex enumeration of baseline data
    double **data_BaseLine = NULL; // baselien data, contains baseline vector, standard deviation

    /*
     * size vertex_enum = data_row x 2
     * size data_BaseLine = data_row x data_column
     */
    vertex_enum = (int **)malloc(data_row * sizeof(int *));
    for (int index = 0; index < data_row; index++)
    {
        *(vertex_enum + index) = (int *)malloc(2 * sizeof(int));
    }
    data_BaseLine = (double **)malloc(data_row * sizeof(double *));
    for (int index = 0; index < data_row; index++)
    {
        *(data_BaseLine + index) = (double *)malloc(data_column * sizeof(double));
    }

    // graph traversal
    for (int index = 0, index_temp = 0; index < var_data->count_vertex && index_temp < data_row; index++)
    {
        AdjListNode *head_node = var_data->array[index].head;

        while (head_node != NULL)
        {
            // vertex enumeration
            vertex_enum[index_temp][0] = index;
            vertex_enum[index_temp][1] = head_node->dest;

            // baseline data
            for (int index_i = 0; index_i < data_column; index_i++)
            {
                data_BaseLine[index_temp][index_i] = head_node->weight[index_i];
            }

            head_node = head_node->next;
            index_temp++;
        }
    }

#if 0 // print data
    puts(">>>>>>>>>>>>vertex_enum:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 2; index_j++)
        {
            printf("%d ", vertex_enum[index_i][index_j]);
        }
        putchar('\n');
    }
    puts("\n>>>>>>>>>>>>data_BaseLine:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            printf("%.4lf ", data_BaseLine[index_i][index_j]);
        }
        for( int index_j = 3; index_j < data_column; index_j++ )
        {
            printf( "%.4le ", data_BaseLine[ index_i ][ index_j ] );
        }
        putchar('\n');
    }
#endif

    // least square equation
    /*
     * sol = \arg\min || Ax - b ||
     * size A = data_row * 3 x cnt_RoverStation * 3
     * size sol = size x = cnt_RoverStation * 3 x 1
     * size b = data_row * 3 x 1
     * residual = Ax - b
     */
    double *lse_sol = NULL, *lse_b = NULL, *lse_residual = NULL;
    double **lse_A = NULL;
    int lse_size_column = var_conf->cnt_RoverStation * 3;
    int lse_size_row = data_row * 3;

    lse_sol = (double *)malloc(lse_size_column * sizeof(double));
    lse_b = (double *)malloc(lse_size_row * sizeof(double));
    lse_residual = (double *)malloc(lse_size_row * sizeof(double));
    lse_A = (double **)malloc(lse_size_row * sizeof(double *));
    for (int index = 0; index < lse_size_row; index++)
    {
        *(lse_A + index) = (double *)malloc(lse_size_column * sizeof(double));
    }

    InitLSE(lse_sol, lse_b, lse_residual, lse_A, lse_size_row, lse_size_column);

    LSESolver_1(var_conf, vertex_enum, data_BaseLine, data_row, data_column, 2,
                lse_A, lse_b, lse_size_row, lse_size_column, lse_sol, lse_residual);

    // updating solution
    LSEUpdateSolution(var_conf, lse_sol, lse_b, lse_size_row, lse_size_column, var_solution);

    // free memory
    for (int index = 0; index < lse_size_row; index++)
    {
        free(*(lse_A + index));
    }
    free(lse_A);
    free(lse_sol);
    free(lse_b);
    free(lse_residual);
    for (int index = 0; index < data_row; index++)
    {
        free(*(data_BaseLine + index));
    }
    free(data_BaseLine);
    for (int index = 0; index < data_row; index++)
    {
        free(*(vertex_enum + index));
    }
    free(vertex_enum);
}

void ImplNetworkAdjustment_2(Config *var_conf, AdjGraph *var_data, Solution *var_solution)
{
    // data size
    int data_row = 0, data_column = 0;
    data_row = var_data->count_edge;
    data_column = sizeof(var_data->array->head->weight) / sizeof(var_data->array->head->weight[0]);

    int **vertex_enum = NULL;      // vertex enumeration of baseline data
    double **data_BaseLine = NULL; // baselien data, contains baseline vector, standard deviation

    /*
     * size vertex_enum = data_row x 2
     * size data_BaseLine = data_row x data_column
     */
    vertex_enum = (int **)malloc(data_row * sizeof(int *));
    for (int index = 0; index < data_row; index++)
    {
        *(vertex_enum + index) = (int *)malloc(2 * sizeof(int));
    }
    data_BaseLine = (double **)malloc(data_row * sizeof(double *));
    for (int index = 0; index < data_row; index++)
    {
        *(data_BaseLine + index) = (double *)malloc(data_column * sizeof(double));
    }

    // graph traversal
    for (int index = 0, index_temp = 0; index < var_data->count_vertex && index_temp < data_row; index++)
    {
        AdjListNode *head_node = var_data->array[index].head;

        while (head_node != NULL)
        {
            // vertex enumeration
            vertex_enum[index_temp][0] = index;
            vertex_enum[index_temp][1] = head_node->dest;

            // baseline data
            for (int index_i = 0; index_i < data_column; index_i++)
            {
                data_BaseLine[index_temp][index_i] = head_node->weight[index_i];
            }

            head_node = head_node->next;
            index_temp++;
        }
    }

#if 0 // print data
    puts(">>>>>>>>>>>>vertex_enum:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 2; index_j++)
        {
            printf("%d ", vertex_enum[index_i][index_j]);
        }
        putchar('\n');
    }
    puts("\n>>>>>>>>>>>>data_BaseLine:");
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            printf("%.4lf ", data_BaseLine[index_i][index_j]);
        }
        for( int index_j = 3; index_j < data_column; index_j++ )
        {
            printf( "%.4le ", data_BaseLine[ index_i ][ index_j ] );
        }
        putchar('\n');
    }
#endif

    // least square equation
    /*
     * sol = \arg\min || Ax - b ||
     * size A = data_row * 3 x cnt_RoverStation * 3
     * size sol = size x = cnt_RoverStation * 3 x 1
     * size b = data_row * 3 x 1
     * residual = Ax - b
     */
    double *lse_sol = NULL, *lse_b = NULL, *lse_residual = NULL;
    double **lse_A = NULL;
    int lse_size_column = var_conf->cnt_RoverStation * 3;
    int lse_size_row = data_row * 3;

    lse_sol = (double *)malloc(lse_size_column * sizeof(double));
    lse_b = (double *)malloc(lse_size_row * sizeof(double));
    lse_residual = (double *)malloc(lse_size_row * sizeof(double));
    lse_A = (double **)malloc(lse_size_row * sizeof(double *));
    for (int index = 0; index < lse_size_row; index++)
    {
        *(lse_A + index) = (double *)malloc(lse_size_column * sizeof(double));
    }

    InitLSE(lse_sol, lse_b, lse_residual, lse_A, lse_size_row, lse_size_column);

    LSESolver_2(var_conf, vertex_enum, data_BaseLine, data_row, data_column, 2,
                lse_A, lse_b, lse_size_row, lse_size_column, lse_sol, lse_residual);

    // updating solution
    LSEUpdateSolution(var_conf, lse_sol, lse_b, lse_size_row, lse_size_column, var_solution);

    // free memory
    for (int index = 0; index < lse_size_row; index++)
    {
        free(*(lse_A + index));
    }
    free(lse_A);
    free(lse_sol);
    free(lse_b);
    free(lse_residual);
    for (int index = 0; index < data_row; index++)
    {
        free(*(data_BaseLine + index));
    }
    free(data_BaseLine);
    for (int index = 0; index < data_row; index++)
    {
        free(*(vertex_enum + index));
    }
    free(vertex_enum);
}

void InitLSE(double *sol, double *rhs, double *residual, double **mat, int row, int column)
{
    // size
    /*
     * residual = mat x sol - rhs
     * size mat = row x column
     * size sol = column x 1
     * size rhs = row x 1
     * size residual = row x 1
     */
    for (int index = 0; index < column; index++)
    {
        *(sol + index) = 0;
    }
    for (int index = 0; index < row; index++)
    {
        *(rhs + index) = 0;
        *(residual + index) = 0;
    }
    for (int index_i = 0; index_i < row; index_i++)
    {
        for (int index_j = 0; index_j < column; index_j++)
        {
            *(*(mat + index_i) + index_j) = 0;
        }
    }
}
