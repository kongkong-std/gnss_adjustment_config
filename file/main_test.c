#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
    FILE *fp = NULL;
    fp = fopen("adjust_data.txt", "rb");

    int node[14][2];
    double data[14][9];

    for (int index = 0; index < 14; index++)
    {
        for (int index_i = 0; index_i < 2; index_i++)
        {
            fscanf(fp, "%d", &node[index][index_i]);
        }
        for (int index_i = 0; index_i < 9; index_i++)
        {
            fscanf(fp, "%lg", &data[index][index_i]);
        }
    }

    fclose(fp);

    // resort
    /*
     * sigma_xx sigma_xy sigma_xz sigma_yy sigma_yz sigma_zz
     *
     * data[][3] sigma_xx
     * data[][4] sigma_yy
     * data[][5] sigma_zz
     * data[][6] sigma_xy
     * data[][7] sigma_yz
     * data[][8] sigma_zx
     */
    double tmp = 0;
    for (int index = 0; index < 14; index++)
    {
        // data[][6] <-> data[][4]
        tmp = data[index][6];
        data[index][6] = data[index][4];
        data[index][4] = tmp;

        // data[][8] <-> data[][5]
        tmp = data[index][8];
        data[index][8] = data[index][5];
        data[index][5] = tmp;
    }

    fp = fopen("adjust_data_sort.txt", "wb");
    for (int index = 0; index < 14; index++)
    {
        for (int index_i = 0; index_i < 2; index_i++)
        {
            fprintf(fp, "%d ", node[index][index_i]);
        }
        for (int index_i = 0; index_i < 3; index_i++)
        {
            fprintf(fp, "%.4lf ", data[index][index_i]);
        }
        for (int index_i = 3; index_i < 9; index_i++)
        {
            fprintf(fp, "%.4le ", data[index][index_i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

#if 0
    for (int index = 0; index < 14; index++)
    {
        for (int index_i = 0; index_i < 2; index_i++)
        {
            printf("%d ", node[index][index_i]);
        }
        for (int index_i = 0; index_i < 3; index_i++)
        {
            printf("%.4lf ", data[index][index_i]);
        }
        for (int index_i = 3; index_i < 9; index_i++)
        {
            printf("%.4le ", data[index][index_i]);
        }
        putchar('\n');
    }
#endif

    return 0;
}

#if 0
int main()
{
    char * buffer = "coo_BaseStation = { 402.35087, -4652995.30109, 4349760.77753, 8086.03178, -4642712.84739, 4360439.08326 }";
    printf( "%s\n", buffer );

    double  * * a = NULL;
    a = ( double * * ) malloc( 2 * sizeof( double * ) );
    for( int index = 0; index < 2; index++ )
    {
	*( a + index ) = ( double * ) malloc( 3 * sizeof( double ) );
    }

#if 0
    for( int index = 0; index < strlen( buffer ); index++ )
    {
	char ch = 0;
	sscanf( buffer + index, "%c", &ch );
	printf( "%c", ch );
    }
#endif

    sscanf( buffer, "%*s = { %lf, %lf, %lf, %lf, %lf, %lf }", *( a ), *( a ) + 1, *( a ) + 2,
	 *( a + 1 ), *( a + 1 ) + 1, *( a + 1) + 2 );

    for( int index_i = 0; index_i < 2; index_i++ )
    {
	for( int index_j = 0; index_j < 3; index_j++ )
	{
	    printf( "%.6lf ", a[ index_i ][ index_j ] );
	}
	putchar( '\n' );
    }

    // free memory
    for( int index = 0; index < 2; index++ )
    {
	free( *( a + index ) );
    }
    free( a );

    return  0;
}
#endif