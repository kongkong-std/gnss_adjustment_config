#include "../include/lse_solver.h"

void update_linsys_full_weight(double **mat_Cov, int row_start, int column_start,
                               double **lse_A, double *lse_b)
{
    double **tmp_mat = NULL;
    double *tmp_vector = NULL;
    tmp_vector = (double *)malloc(3 * sizeof(double));
    tmp_mat = (double **)malloc(3 * sizeof(double *));
    for (int index = 0; index < 3; index++)
    {
        *(tmp_mat + index) = (double *)malloc(3 * sizeof(double));
    }

    // computing tmp_vector = INV(mat_Cov) lse_b( row_start, row_start + 2 )
    for (int index_i = 0; index_i < 3; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            tmp_mat[index_i][index_j] = mat_Cov[index_i][index_j];
        }
    }
    LSEGivens(tmp_mat, lse_b + row_start, 3, 3);
    LSEGauss(tmp_mat, lse_b + row_start, 3, tmp_vector);

    // updating lse_b = tmp_vector
    for (int index = 0; index < 3; index++)
    {
        lse_b[row_start + index] = tmp_vector[index];
    }

    // computing INV( mat_Cov ) lse_A( {row_start, column_start} -- {row_start + 2, column_start + 2} )
    double *tmp_vector_1 = NULL;
    tmp_vector_1 = (double *)malloc(3 * sizeof(double));
    for (int index_j = 0; index_j < 3; index_j++)
    {
        // storing column vector lse_A( row_start : row_start + 2, column_start + index_j )
        for (int index_i = 0; index_i < 3; index_i++)
        {
            tmp_vector_1[index_i] = lse_A[row_start + index_i][column_start + index_j];
        }

        // storing weight matrix
        for (int index_row = 0; index_row < 3; index_row++)
        {
            for (int index_column = 0; index_column < 3; index_column++)
            {
                tmp_mat[index_row][index_column] = mat_Cov[index_row][index_column];
            }
        }

        // solving tmp_mat tmp_vector = tmp_vector_1
        LSEGivens(tmp_mat, tmp_vector_1, 3, 3);
        LSEGauss(tmp_mat, tmp_vector_1, 3, tmp_vector);

        // updating lse_A( row_start : row_start + 2, column_start + index_j ) = tmp_vector
        for (int index_i = 0; index_i < 3; index_i++)
        {
            lse_A[row_start + index_i][column_start + index_j] = tmp_vector[index_i];
        }
    }

    // free memory
    free(tmp_vector_1);
    for (int index = 0; index < 3; index++)
    {
        free(*(tmp_mat + index));
    }
    free(tmp_mat);
    free(tmp_vector);
}

void assemble_mat_Cov(double **data_Baseline, int row_start, double **mat_Cov)
{
    // diagonal element data_BaseLine[][3, 4, 5] sigma_xx, sigma_yy, sigma_zz
    for (int index = 0; index < 3; index++)
    {
        mat_Cov[index][index] = data_Baseline[row_start][3 + index] * data_Baseline[row_start][3 + index];
    }

    // sub-diagonal element data_BaseLine[][6, 7] sigma_xy, sigma_yz
    for (int index = 0; index < 2; index++)
    {
        if (data_Baseline[row_start][6 + index] > 0)
        {
            mat_Cov[index][index + 1] = data_Baseline[row_start][6 + index] * data_Baseline[row_start][6 + index];
            mat_Cov[index + 1][index] = data_Baseline[row_start][6 + index] * data_Baseline[row_start][6 + index];
        }
        else
        {
            mat_Cov[index][index + 1] = -data_Baseline[row_start][6 + index] * data_Baseline[row_start][6 + index];
            mat_Cov[index + 1][index] = -data_Baseline[row_start][6 + index] * data_Baseline[row_start][6 + index];
        }
    }

    // up-right angular element data_BaseLine[][8] sigma_zx
    if (data_Baseline[row_start][8] > 0)
    {
        mat_Cov[0][2] = data_Baseline[row_start][8] * data_Baseline[row_start][8];
        mat_Cov[2][0] = data_Baseline[row_start][8] * data_Baseline[row_start][8];
    }
    else
    {
        mat_Cov[0][2] = -data_Baseline[row_start][8] * data_Baseline[row_start][8];
        mat_Cov[2][0] = -data_Baseline[row_start][8] * data_Baseline[row_start][8];
    }
}

void LSEGauss(double **mat, double *rhs, int column, double *sol)
{
    // solving upper triangular linear system with direct method
    /*
     * x_k = ( b_k - \sum{ j = k + 1, n } a_kj x_j ) / a_kk
     * */
    for (int index = column - 1; index >= 0; index--)
    {
        double sum = 0;
        for (int index_j = index + 1; index_j < column; index_j++)
        {
            sum += mat[index][index_j] * sol[index_j];
        }
        sol[index] = (rhs[index] - sum) / mat[index][index];
    }
}

void givens_rotation_impl(double val_a, double val_b, int index_i, int index_j,
                          int row, int column, double **mat, double *rhs)
{
    double givens_tau = 0, givens_c = 0, givens_s = 0;
    if (val_a != 0 && val_b != 0)
    {
        if (fabs(val_b) > fabs(val_a))
        {
            givens_tau = -val_a / val_b;
            givens_s = 1 / sqrt(1 + givens_tau * givens_tau);
            givens_c = givens_s * givens_tau;
        }
        else
        {
            givens_tau = -val_b / val_a;
            givens_c = 1 / sqrt(1 + givens_tau * givens_tau);
            givens_s = givens_c * givens_tau;
        }
    }
    else if (val_a != 0 && val_b == 0)
    {
        givens_c = 1;
        givens_s = 0;
    }
    else if (val_a == 0 && val_b != 0)
    {
        givens_c = 0;
        givens_s = 1;
    }
    else if (val_a == 0 && val_b == 0)
    {
        givens_c = 1;
        givens_s = 0;
    }

    // computing mat
    /*
     * mat( index_i, index ) = c * mat( index_i, index ) - s * mat( index_j, index )
     * mat( index_j, index ) = s * mat( index_i, index ) + c * mat( index_j, index )
     * */
    double temp_1 = 0, temp_2 = 0;
    for (int index = 0; index < column; index++)
    {
        temp_1 = mat[index_i][index];
        temp_2 = mat[index_j][index];
        mat[index_i][index] = givens_c * temp_1 - givens_s * temp_2;
        mat[index_j][index] = givens_s * temp_1 + givens_c * temp_2;
    }

    // computing rhs
    /*
     * rhs( index_i, 1 ) = c * rhs( index_i, 1 ) - s * rhs( index_j, 1 )
     * rhs( index_j, 1 ) = s * rhs( index_i, 1 ) + c * rhs( index_j, 1 )
     * */
    temp_1 = rhs[index_i];
    temp_2 = rhs[index_j];
    rhs[index_i] = givens_c * temp_1 - givens_s * temp_2;
    rhs[index_j] = givens_s * temp_1 + givens_c * temp_2;
}

void LSEGivens(double **mat, double *rhs, int row, int column)
{
    for (int index = 0; index < column; index++)
    {
        for (int index_j = index + 1; index_j < row; index_j++)
        {
            givens_rotation_impl(mat[index][index], mat[index_j][index], index, index_j,
                                 row, column, mat, rhs);
        }
    }
}

int IsBaseStation(int code, int *code_BaseStation, int cnt_BaseStation)
{
    int value = 0;

    for (int index = 0; index < cnt_BaseStation; index++)
    {
        if (*(code_BaseStation + index) == code)
        {
            value = 1;
            break;
        }
    }

    return value;
}

void LSEAssemble_2(int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                   int *code_BaseStation, double **coo_BaseStation, int cnt_BaseStation, int coo_column,
                   double **lse_A, double *lse_b, int lse_size_row, int lse_size_column)
{
    // full weight W = INV{ Cov }
    /*
     * WAx - Wb = INV( mat_Cov )Ax - INV( mat_Cov ) b
     */
    double **mat_Cov = NULL; // size mat_Cov = 3 x 3
    mat_Cov = (double **)malloc(3 * sizeof(double *));
    for (int index = 0; index < 3; index++)
    {
        *(mat_Cov + index) = (double *)malloc(3 * sizeof(double));
    }

    for (int index_i = 0; index_i < data_row; index_i++)
    {
        if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
            IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // baseline: base station -> base station
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // base line: base station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * P _ rover = P _ base + P _ baseline
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] = 1.;
                lse_b[3 * index_i + index_j] = coo_BaseStation[vertex_enum[index_i][0]][index_j] +
                                               data_BaseLine[index_i][index_j];
            }

            // assemble mat_Cov
            assemble_mat_Cov(data_BaseLine, index_i, mat_Cov);

            // updating linear system
            /*
             * lse_A = INV( mat_Cov ) lse_A
             * lse_b = INV( mat_Cov ) lse_b
             */
            update_linsys_full_weight(mat_Cov, 3 * index_i, (vertex_enum[index_i][1] - cnt_BaseStation) * 3,
                                      lse_A, lse_b);
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // base line: rover station -> base station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * -P_rover = P _ baseline - P _ base
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] = -1.;
                lse_b[3 * index_i + index_j] = data_BaseLine[index_i][index_j] -
                                               coo_BaseStation[vertex_enum[index_i][1]][index_j];
            }

            // assemble mat_Cov
            assemble_mat_Cov(data_BaseLine, index_i, mat_Cov);

            // updating linear system
            /*
             * lse_A = INV( mat_Cov ) lse_A
             * lse_b = INV( mat_Cov ) lse_b
             */
            update_linsys_full_weight(mat_Cov, 3 * index_i, (vertex_enum[index_i][0] - cnt_BaseStation) * 3,
                                      lse_A, lse_b);
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // base line: rover station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side:
             * P _ rover_1 - P _ rover_2 = P _ baseline
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] = 1.;
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] = -1.;
                lse_b[3 * index_i + index_j] = data_BaseLine[index_i][index_j];
            }

            // assemble mat_Cov
            assemble_mat_Cov(data_BaseLine, index_i, mat_Cov);

            /*
             * updating linear system
             * lse_A = INV( mat_Cov ) lse_A
             * lse_b = INV( mat_Cov ) lse_b
             */
            update_linsys_full_weight(mat_Cov, 3 * index_i, (vertex_enum[index_i][1] - cnt_BaseStation) * 3,
                                      lse_A, lse_b);
            update_linsys_full_weight(mat_Cov, 3 * index_i, (vertex_enum[index_i][0] - cnt_BaseStation) * 3,
                                      lse_A, lse_b);
        }
    }

    // free memory
    for (int index = 0; index < 3; index++)
    {
        free(*(mat_Cov + index));
    }
    free(mat_Cov);
}

void LSEAssemble_1(int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                   int *code_BaseStation, double **coo_BaseStation, int cnt_BaseStation, int coo_column,
                   double **lse_A, double *lse_b, int lse_size_row, int lse_size_column)
{
    // diagonal weight W = INV{ diag(covariance) }
    /*
     * mat_Cov = diag(covariance)
     * WAx - Wb = INV( mat_Cov ) Ax - INV( mat_Cov ) b
     */
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
            IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // baseline: base station -> base station
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // base line: base station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * P _ rover = P _ base + P _ baseline
             * lse_A / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             * lse_b / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] =
                    1. / data_BaseLine[index_i][3 + index_j] / data_BaseLine[index_i][3 + index_j];
                lse_b[3 * index_i + index_j] = (coo_BaseStation[vertex_enum[index_i][0]][index_j] +
                                                data_BaseLine[index_i][index_j]) /
                                               data_BaseLine[index_i][3 + index_j] /
                                               data_BaseLine[index_i][3 + index_j];
            }
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // baseline: rover station -> base station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * -P_rover = P _ baseline - P _ base
             * lse_A / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             * lse_b / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] =
                    -1. / data_BaseLine[index_i][3 + index_j] / data_BaseLine[index_i][3 + index_j];
                lse_b[3 * index_i + index_j] = (data_BaseLine[index_i][index_j] -
                                                coo_BaseStation[vertex_enum[index_i][1]][index_j]) /
                                               data_BaseLine[index_i][3 + index_j] /
                                               data_BaseLine[index_i][3 + index_j];
            }
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // baseline: rover station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side:
             * P _ rover_1 - P _ rover_2 = P _ baseline
             * lse_A / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             * lse_b / { data_BaseLine[index_i][3, 4, 5] ^ 2 }
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] =
                    1. / data_BaseLine[index_i][3 + index_j] / data_BaseLine[index_i][3 + index_j];
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] =
                    -1. / data_BaseLine[index_i][3 + index_j] / data_BaseLine[index_i][3 + index_j];
                lse_b[3 * index_i + index_j] = data_BaseLine[index_i][index_j] /
                                               data_BaseLine[index_i][3 + index_j] / data_BaseLine[index_i][3 + index_j];
            }
        }
    }
}

void LSEAssemble_0(int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                   int *code_BaseStation, double **coo_BaseStation, int cnt_BaseStation, int coo_column,
                   double **lse_A, double *lse_b, int lse_size_row, int lse_size_column)
{
    // equal weight W = I
    /*
     * WAx - Wb = Ax - b
     */
    for (int index_i = 0; index_i < data_row; index_i++)
    {
        if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
            IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // baseline: base station -> base station
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 1 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // base line: base station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * P _ rover = P _ base + P _ baseline
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] = 1.;
                lse_b[3 * index_i + index_j] = coo_BaseStation[vertex_enum[index_i][0]][index_j] +
                                               data_BaseLine[index_i][index_j];
            }
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 1)
        {
            // base line: rover station -> base station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side: 3*index_i -- 3*index_i + 2
             * -P_rover = P _ baseline - P _ base
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] = -1.;
                lse_b[3 * index_i + index_j] = data_BaseLine[index_i][index_j] -
                                               coo_BaseStation[vertex_enum[index_i][1]][index_j];
            }
        }
        else if (IsBaseStation(vertex_enum[index_i][0], code_BaseStation, cnt_BaseStation) == 0 &&
                 IsBaseStation(vertex_enum[index_i][1], code_BaseStation, cnt_BaseStation) == 0)
        {
            // base line: rover station -> rover station
            // coefficient matrix and right-hand side vector
            /*
             * { 3*index_i, ([index_i][1] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][1] - cnt_base)*3 + 2 } diag(1)
             * { 3*index_i, ([index_i][0] - cnt_base)*3 } -- { 3*index_i + 2, ([index_i][0] - cnt_base)*3 + 2 } diag(-1)
             * right-hand side:
             * P _ rover_1 - P _ rover_2 = P _ baseline
             */
            for (int index_j = 0; index_j < 3; index_j++)
            {
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][1] - cnt_BaseStation) * 3 + index_j] = 1.;
                lse_A[3 * index_i + index_j][(vertex_enum[index_i][0] - cnt_BaseStation) * 3 + index_j] = -1.;
                lse_b[3 * index_i + index_j] = data_BaseLine[index_i][index_j];
            }
        }
    }
}

void LSESolver_0(Config *var_conf, int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                 double **lse_A, double *lse_b, int lse_size_row, int lse_size_column, double *lse_sol, double *lse_residual)
{
    int *code_BaseStation = NULL; // base station code, size = cnt_BaseStation x 1
    code_BaseStation = (int *)malloc(var_conf->cnt_BaseStation * sizeof(int));

    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(code_BaseStation + index) = index;
    }

    double **coo_BaseStation = NULL; // base station position, size = cnt_BaseStation x 3
    coo_BaseStation = (double **)malloc(var_conf->cnt_BaseStation * sizeof(double *));
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(coo_BaseStation + index) = (double *)malloc(3 * sizeof(double));
    }
    for (int index_i = 0; index_i < var_conf->cnt_BaseStation; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            *(*(coo_BaseStation + index_i) + index_j) = var_conf->coo_BaseStation[3 * index_i + index_j];
        }
    }

    // assemble least square equation
    LSEAssemble_0(vertex_enum, data_BaseLine, data_row, data_column, vertex_column,
                  code_BaseStation, coo_BaseStation, var_conf->cnt_BaseStation, 3,
                  lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // givens rotation transform least square equation
    LSEGivens(lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE after givens rotation
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // solving uptriangular linear system
    LSEGauss(lse_A, lse_b, lse_size_column, lse_sol);

#if 0 // display solution
    puts(">>>>>>>>>>>> lse_sol:");
    for (int index = 0; index < lse_size_column; index++)
    {
        printf("%.6lf\n", lse_sol[index]);
    }
#endif

    // free memory
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        free(*(coo_BaseStation + index));
    }
    free(coo_BaseStation);
    free(code_BaseStation);
}

void LSESolver_1(Config *var_conf, int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                 double **lse_A, double *lse_b, int lse_size_row, int lse_size_column, double *lse_sol, double *lse_residual)
{
    int *code_BaseStation = NULL; // base station code, size = cnt_BaseStation x 1
    code_BaseStation = (int *)malloc(var_conf->cnt_BaseStation * sizeof(int));

    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(code_BaseStation + index) = index;
    }

    double **coo_BaseStation = NULL; // base station position, size = cnt_BaseStation x 3
    coo_BaseStation = (double **)malloc(var_conf->cnt_BaseStation * sizeof(double *));
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(coo_BaseStation + index) = (double *)malloc(3 * sizeof(double));
    }
    for (int index_i = 0; index_i < var_conf->cnt_BaseStation; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            *(*(coo_BaseStation + index_i) + index_j) = var_conf->coo_BaseStation[3 * index_i + index_j];
        }
    }

    // assemble least square equation
    LSEAssemble_1(vertex_enum, data_BaseLine, data_row, data_column, vertex_column,
                  code_BaseStation, coo_BaseStation, var_conf->cnt_BaseStation, 3,
                  lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // givens rotation transform least square equation
    LSEGivens(lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE after givens rotation
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // solving uptriangular linear system
    LSEGauss(lse_A, lse_b, lse_size_column, lse_sol);

#if 0 // display solution
    puts(">>>>>>>>>>>> lse_sol:");
    for (int index = 0; index < lse_size_column; index++)
    {
        printf("%.6lf\n", lse_sol[index]);
    }
#endif

    // free memory
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        free(*(coo_BaseStation + index));
    }
    free(coo_BaseStation);
    free(code_BaseStation);
}

void LSESolver_2(Config *var_conf, int **vertex_enum, double **data_BaseLine, int data_row, int data_column, int vertex_column,
                 double **lse_A, double *lse_b, int lse_size_row, int lse_size_column, double *lse_sol, double *lse_residual)
{
    int *code_BaseStation = NULL; // base station code, size = cnt_BaseStation x 1
    code_BaseStation = (int *)malloc(var_conf->cnt_BaseStation * sizeof(int));

    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(code_BaseStation + index) = index;
    }

    double **coo_BaseStation = NULL; // base station position, size = cnt_BaseStation x 3
    coo_BaseStation = (double **)malloc(var_conf->cnt_BaseStation * sizeof(double *));
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        *(coo_BaseStation + index) = (double *)malloc(3 * sizeof(double));
    }
    for (int index_i = 0; index_i < var_conf->cnt_BaseStation; index_i++)
    {
        for (int index_j = 0; index_j < 3; index_j++)
        {
            *(*(coo_BaseStation + index_i) + index_j) = var_conf->coo_BaseStation[3 * index_i + index_j];
        }
    }

    // assemble least square equation
    LSEAssemble_2(vertex_enum, data_BaseLine, data_row, data_column, vertex_column,
                  code_BaseStation, coo_BaseStation, var_conf->cnt_BaseStation, 3,
                  lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // givens rotation transform least square equation
    LSEGivens(lse_A, lse_b, lse_size_row, lse_size_column);

#if 0 // display LSE after givens rotation
    printf(">>>>>>>>>>>> lse_size_row = %d, lse_size_column = %d\n", lse_size_row, lse_size_column);
    puts(">>>>>>>>>>>> lse_A:");
    for (int index_i = 0; index_i < lse_size_row; index_i++)
    {
        for (int index_j = 0; index_j < lse_size_column; index_j++)
        {
            printf("%.1lf ", lse_A[index_i][index_j]);
        }
        putchar('\n');
    }
    puts(">>>>>>>>>>>> lse_b:");
    for (int index = 0; index < lse_size_row; index++)
    {
        printf("%.4lf\n", lse_b[index]);
    }
#endif

    // solving uptriangular linear system
    LSEGauss(lse_A, lse_b, lse_size_column, lse_sol);

#if 0 // display solution
    puts(">>>>>>>>>>>> lse_sol:");
    for (int index = 0; index < lse_size_column; index++)
    {
        printf("%.6lf\n", lse_sol[index]);
    }
#endif

    // free memory
    for (int index = 0; index < var_conf->cnt_BaseStation; index++)
    {
        free(*(coo_BaseStation + index));
    }
    free(coo_BaseStation);
    free(code_BaseStation);
}
