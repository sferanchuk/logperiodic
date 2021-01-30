
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double linear_regression( double *x, double *y, int dim, double *slope, double *intercept )
{
    int i;
    double sumx = 0, sumy = 0, sumxx = 0, sumxy = 0, sumyy = 0;
    
    for ( i = 0; i < dim; i++ )
    {
        sumx += x[i];
        sumy += y[i];
        sumxx += x[i] * x[i];
        sumxy += x[i] * y[i];
        sumyy += y[i] * y[i];
    }
    *slope = ( dim * sumxy - sumx * sumy ) / ( dim * sumxx - sumx * sumx );
    *intercept = sumy / dim - sumx * (*slope) / dim;
    return ( dim * sumxy - sumx * sumy ) / sqrt( ( dim * sumxx - sumx * sumx ) * ( dim * sumyy - sumy * sumy ) );
}

double path_length( double *data, double dx, int dim )
{
    double sum = 0, diff;
    int i;
    double dx2 = dx * dx;
    
    for ( i = 0; i + 1 < dim; i++ )
    {
        diff = data[ i + 1 ] - data[i];
        sum += sqrt( diff * diff + dx2 );
    }
    return sum;
}

void higuchi_dimension( double *data, int npoints, double *slope, double *intercept, double *correlation )
{
    double *log_path, *log_interval;
    int i, j, k;
    int interval, slice_dim, num_logpoints;
    double *data_slice;
    double partial_sum;
    
    num_logpoints = (int) log( npoints );
    log_path = (double*) malloc( num_logpoints * sizeof(double) );
    log_interval = (double*) malloc( num_logpoints * sizeof( double ) );
    
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( i );
        slice_dim = npoints / interval;
        data_slice = (double*) malloc( slice_dim * sizeof(double) );
        partial_sum = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; k < slice_dim; k++ )
            {
                data_slice[k] = data[ j + interval * k ];
            }
            partial_sum += path_length( data_slice, 0, slice_dim );
        }
        log_path[i] = -log( ( partial_sum / interval ) * ( (double)npoints / ( interval * slice_dim ) ) );
        log_interval[i] = i;
        free( data_slice );
    }
    *correlation = linear_regression( log_interval, log_path, num_logpoints, slope, intercept );
    free( log_interval );
    free( log_path );
}

/*
double path_length_xy( double *datax, double *datay, int dim )
{
    double sum = 0, diffx, diffy;
    int i;
    
    for ( i = 0; i + 1 < dim; i++ )
    {
        diffx = datax[ i + 1 ] - datax[i];
        diffy = datay[ i + 1 ] - datay[i];
        sum += sqrt( diffx * diffx + diffy * diffy );
    }
    return sum;
}

double path_dimension( double *data, int npoints )
{
    double *log_path, *log_interval;
    int i, j, k;
    int interval, slice_dim, num_logpoints;
    double *data_slice, data_min, data_max;
    double partial_sum, fdim_slope, fdim_intercept;
    
    data_min = data_max = data[0]; for ( i = 1; i < npoints; i++ ) { if ( data[i] > data_max ) data_max = data[i]; if ( data[i] < data_min ) data_min = data[i]; } 

    num_logpoints = (int) log( npoints );
    log_path = (double*) malloc( num_logpoints * sizeof(double) );
    log_interval = (double*) malloc( num_logpoints * sizeof( double ) );
    
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( i );
        slice_dim = npoints / interval;
        data_slice = (double*) malloc( slice_dim * sizeof(double) );
        partial_sum = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; k < slice_dim; k++ )
            {
                data_slice[k] = data[ j + interval * k ];
            }
            partial_sum += path_length( data_slice, interval * ( data_max - data_min ) / npoints, slice_dim );
        }
        log_path[i] = log( ( partial_sum / interval ) * ( (double)npoints / ( interval * slice_dim ) ) );
        log_interval[i] = i;
        free( data_slice );
    }
    linear_regression( log_interval, log_path, num_logpoints, &fdim_slope, &fdim_intercept );
    free( log_interval );
    free( log_path );
    return -fdim_slope;
}

double path_dimension_xy( double *data_x, double *data_y, int npoints )
{
    double *log_path, *log_interval;
    int i, j, k;
    int interval, slice_dim, num_logpoints;
    double *slice_x, *slice_y, y_max, y_min, scale;
    double partial_sum, fdim_slope, fdim_intercept;
    
    y_min = y_max = data_y[0]; for ( i = 1; i < npoints; i++ ) { if ( data_y[i] > y_max ) y_max = data_y[i]; if ( data_y[i] < y_min ) y_min = data_y[i]; } 
    scale = ( y_max - y_min ) / ( data_x[ npoints - 1 ] - data_x[ 0 ] );
    num_logpoints = (int) log( npoints );
    log_path = (double*) malloc( num_logpoints * sizeof(double) );
    log_interval = (double*) malloc( num_logpoints * sizeof( double ) );
   
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( i );
        slice_dim = npoints / interval;
        slice_x = (double*) malloc( slice_dim * sizeof(double) );
        slice_y = (double*) malloc( slice_dim * sizeof(double) );
        partial_sum = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; k < slice_dim; k++ )
            {
                slice_x[k] = data_x[ j + interval * k ] * scale;
                slice_y[k] = data_y[ j + interval * k ];
            }
            partial_sum += path_length_xy( slice_x, slice_y, slice_dim );
        }
        log_path[i] = log( ( partial_sum / interval ) * ( (double)npoints / ( interval * slice_dim ) ) );
        log_interval[i] = i;
        free( slice_x );
        free( slice_y );
    }
    linear_regression( log_interval, log_path, num_logpoints, &fdim_slope, &fdim_intercept );
    free( log_interval );
    free( log_path );
    return -fdim_slope;
}
*/

void peng_dimension( double *data_x, double *data_y, int npoints, double *slope, double *intercept, double *correlation )
{
    double *log_path, *log_interval;
    int i, j, k;
    int interval, num_slices, num_logpoints;
    double *slice_x, *slice_y;
    double slice_slope, slice_intercept, sum_slices;
    
    num_logpoints = (int) log( npoints );
    log_path = (double*) malloc( num_logpoints * sizeof(double) );
    log_interval = (double*) malloc( num_logpoints * sizeof( double ) );
    
    for ( i = 2; i < num_logpoints; i++ )
    {
        interval = (int) exp( i );
        num_slices = npoints / interval;
        slice_x = (double*) malloc( interval * sizeof(double) );
        slice_y = (double*) malloc( interval * sizeof(double) );
        sum_slices = 0;
        for ( j = 0; j < num_slices; j++ )
        {
            for ( k = 0; k < interval; k++ )
            {
                slice_x[k] = data_x[ j * interval + k ];
                slice_y[k] = data_y[ j * interval + k ];
            }
            linear_regression( slice_x, slice_y, interval, &slice_slope, &slice_intercept );
            sum_slices += ( slice_x[ interval - 1 ] - slice_x[0] ) * ( 1 + slice_slope * slice_slope );
        }
        log_path[i] = log( ( sum_slices / num_slices ) * ( (double)npoints / ( interval * num_slices ) ) );
        log_interval[i] = i;
        free( slice_x );
        free( slice_y );
    }
    *correlation = linear_regression( log_interval, log_path, num_logpoints, slope, intercept );
    free( log_interval );
    free( log_path );
}

#define MAX_OUTBREAKS 10
#define MAX_PERIODS 100
#define MAX_SAMPLES 1000

static double abs_diff( double p1, double p2 ) { if ( p1 == p2 ) return 0; return ( p1 - p2 > 0 ) ? ( p1 - p2 ) : ( p2 - p1 ); } 
static int straight_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? 1 : -1; } 
static int reversed_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? -1 : 1; } 

void logperiodic_approximation( double *data_x, double *data_y, int num_points, int print_datapoints, int *res_direction, int *res_period,  double *res_slope, double *res_intercept, double *res_correlation )
{
    int i, j, k, num_samples, num_outbreaks, num_periods;
    double sum, sum_sq, d, vmin, slope, intercept, eps;
    double global_mean, local_mean, global_sd, local_sd, outbreak_bound;
    double log_distr[ MAX_SAMPLES ], log_count[ MAX_SAMPLES ];
    int indices[ MAX_OUTBREAKS] ;
    double outbreak_x[ MAX_OUTBREAKS ], outbreak_y[ MAX_OUTBREAKS ], *ptr_x, *ptr_y;
    double periods[2][ MAX_PERIODS ];

    num_samples = ( num_points < MAX_SAMPLES ) ? num_points : MAX_SAMPLES;
    sum = 0;
    sum_sq = 0;
    for ( i = 0; i < num_samples; i++ )
    {
        sum += data_y[i];
        sum_sq += data_y[i] * data_y[i];
    }
    global_mean = sum / num_samples;
    global_sd = sqrt( sum_sq / num_samples - global_mean * global_mean );
    eps = 0.0001 * global_sd;
    for ( i = 0; i < num_samples; i++ ) log_count[i] = log( i + 1 );
    for ( i = 0; i < num_samples; i++ ) log_distr[i] = log( ( abs_diff( data_y[i], global_mean ) < eps ) ? eps : abs_diff( data_y[i], global_mean ) );
    qsort( log_distr, i, sizeof( double), reversed_compare );
    linear_regression( log_count, log_distr, i, &slope, &intercept );
    outbreak_bound = exp( fabs( intercept ) );
 
    for ( ptr_y = outbreak_y, ptr_x = outbreak_x, num_outbreaks = 0; num_outbreaks < MAX_OUTBREAKS && *ptr_y < outbreak_bound; ptr_y++, ptr_x++, num_outbreaks++ )
    {
        *ptr_y = 0;
        for ( i = 0; i < num_points; i++ )
        {
            for ( j = 0; j < num_outbreaks; j++ ) if ( indices[j] == i ) break; 
            if ( j < num_outbreaks ) continue;
            
            vmin = abs_diff( data_y[i], global_mean ); 
            for ( j = 3; j + i < num_points && j < i; j *= 2 )
            {
                sum = 0;
                sum_sq = 0;
                for ( k = i - j; k <= i + j; k++ )
                {
                    sum += data_y[ k ];
                    sum_sq += data_y[k] * data_y[k];
                }
                local_mean = sum / ( 2 * j + 1 );
                local_sd = sqrt( sum_sq / ( 2 * j + 1 ) - local_mean * local_mean );
                d = abs_diff( data_y[i], local_mean ) * global_sd / local_sd;
                if ( d < vmin ) vmin = d;
            }
            if ( vmin >= *ptr_y ) 
            {
                *ptr_y = vmin;
                *ptr_x = data_x[i];
                indices[ num_outbreaks ] = i;
            }
        }
    }
    
    qsort( outbreak_x, num_outbreaks, sizeof( double ), straight_compare );

    num_periods = MAX_PERIODS - ( MAX_OUTBREAKS - num_outbreaks );
    periods[0][0] = periods[1][0] = 0; 
    for ( i = 1; i < num_periods; i++ ) periods[0][i] = log( i + 1 ) + periods[0][ i - 1 ];
    for ( i = 1; i < num_periods; i++ ) periods[1][i] = log( num_periods - i ) + periods[1][ i - 1 ];

    vmin = -1;
    for ( i = 0; i < num_periods - num_outbreaks; i++ )
    {
        for ( j = 0; j < 2; j++ )
        {
            d = linear_regression( periods[j] + i, outbreak_x, num_outbreaks, &slope, &intercept );
            if ( vmin < 0 || vmin > d )
            {
                vmin = d;
                *res_direction = j;
                *res_period = i;   
            }
        }
    }
    
    *res_correlation = linear_regression( periods[ *res_direction ] + *res_period, outbreak_x, num_outbreaks, res_slope, res_intercept );
 
    if ( print_datapoints ) 
    {
        printf( "outbreak_x: " );
        for ( j = 0; j < num_outbreaks; j++ ) printf( "%7g ", outbreak_x[j] );
        printf( "\n" );
        printf( "outbreak_y: " );
        for ( i = 0; i < num_outbreaks; i++ ) 
        for ( j = 0; j < num_outbreaks; j++ )
        {
            if ( data_x[ indices[j] ] == outbreak_x[i] ) printf( "%7g ", data_y[ indices[j] ] ); 
        }
        printf( "\n" );
        printf( "approximated periods: " );
        for ( j = 0; j < num_outbreaks; j++ ) printf( "%7g ", periods[ *res_direction ][ *res_period  + j ] * ( *res_slope ) + *( res_intercept ) );
        printf( "\n" );
    }
}

int main( int argc, char **argv )
{
    FILE *ifile = stdin;
    char buf[1024] = "";
    char *p;
    int i, maxnpoints = 0;
    int npoints = 0;
    double *values = 0, *ptr, *xvalues, yvalues;
    int print_outbreaks;
    double res_correlation, res_dimension, res_intercept;
    int res_direction, res_period;
    
    if ( argc == 1 ) 
    { 
        printf( "arguments: -h (higuchi fd) | -p (peng f.d.) | -x (peng_xy) | -c (log-periodic);\n -d to print outbreaks in log-periodc model;\nstdin: input data, in one line\n" ); return 1;
    }
    
    p=buf;
    while ( fread( p, 1, 1, ifile ) == 1 )
    {
        if ( *p == ',' )
        {
            if ( npoints >= maxnpoints - 1 )
            {
                maxnpoints += 50;
                ptr = (double*) malloc( maxnpoints * sizeof( double ) );
                memcpy( ptr, values, npoints * sizeof( double ) );
                if ( values != 0 ) free( values );
                values = ptr;
            }
            *(p + 1) = 0;
            values[ npoints++ ] = atof( buf );
            p = buf;
        }
        else 
        {
            p++;
            if ( p - buf >= 1024 ) { printf( "buffer overflow\n" ); return 1; }
        }
    }
    values[ npoints++ ] = atof( buf );
    if ( npoints > 0 && argv[1][0] == '-' ) 
    {
        switch ( argv[1][1] )
        {
            case 'h':
            {
                higuchi_dimension( values, npoints, &res_dimension, &res_intercept, &res_correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", res_dimension, res_intercept, res_correlation );        
                free( values );
                break;
            }
            case 'p':
            {
                double *xvalues = (double*) malloc( sizeof( double ) * npoints );
                for ( i = 0; i < npoints; i++ ) xvalues[i] = i;
                peng_dimension( xvalues, values, npoints,  &res_dimension, &res_intercept, &res_correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", res_dimension, res_intercept, res_correlation );        
                free( values );
                free( xvalues );
                break;
            }        
            case 'x':
            {
                npoints /= 2;
                double *xvalues = (double*) malloc( sizeof( double ) * npoints );
                double *yvalues = (double*) malloc( sizeof( double ) * npoints );
                for ( i = 0; i < npoints; i++ )
                {
                    xvalues[i] = values[ 2 * i ];
                    yvalues[i] = values[ 2 * i + 1 ];
                }
                peng_dimension( xvalues, yvalues, npoints, &res_dimension, &res_intercept, &res_correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", res_dimension, res_intercept, res_correlation );        
                free( values );
                free( xvalues );
                free( yvalues );
                break;
            }
            case 'c':
            case 'd':
            {
                npoints /= 2;
                double *xvalues = (double*) malloc( sizeof( double ) * npoints );
                double *yvalues = (double*) malloc( sizeof( double ) * npoints );
                for ( i = 0; i < npoints; i++ )
                {
                    xvalues[i] = values[ 2 * i ];
                    yvalues[i] = values[ 2 * i + 1 ];
                }
                print_outbreaks = ( argv[1][1] == 'c' ) ? 0 : 1;
                logperiodic_approximation( xvalues, yvalues, npoints, print_outbreaks, &res_direction, &res_period, &res_dimension, &res_intercept, &res_correlation );
                printf( "direction: %s\nperiod in beginning: %d\nrate: %g\nphase: %g\ncorrelation: %g\n", ( res_direction == 0 ) ? "expansion from past" : "contraction to future", res_period, res_dimension, res_intercept, res_correlation );          
                free( values );
                free( xvalues );
                free( yvalues );
                break;
            }
            default:
                printf( "incorrect option key %s\n", argv[1] );
                free( values );
                break;
        }
    }
    return 0;
}
