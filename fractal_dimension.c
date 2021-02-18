
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

#define MAX_INTERVALS 1000

void higuchi_dimension( double *data_x, double *data_y, int npoints, int method, int interval_min, int interval_max, int extra_output, int linincr, double *slope, double *intercept, double *correlation )
{
    int i, j, k, npasses;
    int num_intervals, interval;
    double dx, dy, partial_sum;
    double path[ MAX_INTERVALS ], intervals[ MAX_INTERVALS ], log_path[ MAX_INTERVALS ], log_intervals[ MAX_INTERVALS ];
    
    if ( linincr ) 
    { 
       for ( i = ( ( interval_min == 0 ) ? 2 : interval_min ), num_intervals = 0; i < npoints && num_intervals < MAX_INTERVALS && ( interval_max == 0 || i < interval_max ); i++, num_intervals++ ) intervals[ num_intervals ] = i; 
    }
    else 
    { 
        for ( i = ( ( interval_min == 0 ) ? 2 : interval_min ), num_intervals = 0; i < npoints && ( interval_max == 0 || i < interval_max ); i = ( 3 * i ) / 2, num_intervals++ ) intervals[ num_intervals] = i; 
    }

    for ( i = 0; i < num_intervals; i++ )
    {
        interval = intervals[i];
        partial_sum = 0;
        npasses = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; j + interval * ( k + 1 ) < npoints; k++ )
            {
                dx = data_x[ interval * ( k + 1 ) + j ] - data_x[ interval * k + j ];
                dy = data_y[ interval * ( k + 1 ) + j ] - data_y[ interval * k + j ];
                partial_sum += ( method == 0 ) ? ( ( dy > 0 ) ? dy : -dy ) : sqrt( dx * dx + dy * dy );
                npasses++;
            }
        }
        path[i] = ( partial_sum / npasses ) * ( (double)npoints / interval );
    }
    for ( i = 0; i < num_intervals; i++ )
    {
            log_intervals[i] = log( intervals[i] );
            log_path[i] = log( path[i] );
    }
    *correlation = linear_regression( log_intervals, log_path, num_intervals, slope, intercept );
    
    if ( extra_output )
    {
        printf( "path:" );
        for ( i = 0; i < num_intervals; i++ ) printf( " %g,%g", log_intervals[i], log_path[i] );
        printf( "\n" );
    }
}

void peng_dimension( double *data_x, double *data_y, int npoints, int method, int interval_min, int interval_max, int extra_output, int linincr, double *slope, double *intercept, double *correlation )
{
    int i, j, k, npasses;
    int num_intervals, interval;
    double slice_slope, slice_intercept, dx, dy, partial_sum;
    double path[ MAX_INTERVALS ], intervals[ MAX_INTERVALS ], log_path[ MAX_INTERVALS ], log_intervals[ MAX_INTERVALS ];
  
    if ( linincr ) 
    { 
       for ( i = ( ( interval_min == 0 ) ? 3 : interval_min ), num_intervals = 0; i < npoints && num_intervals < MAX_INTERVALS && ( interval_max == 0 || i < interval_max ); i++, num_intervals++ ) intervals[ num_intervals ] = i; 
    }
    else 
    { 
        for ( i = ( ( interval_min == 0 ) ? 3 : interval_min ), num_intervals = 0; i < npoints && ( interval_max == 0 || i < interval_max ); i = ( 3 * i ) / 2, num_intervals++ ) intervals[ num_intervals] = i; 
    }
    
    for ( i = 0; i < num_intervals; i++ )
    {
        interval = intervals[i];
        partial_sum = 0;
        npasses = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; j + interval * ( k + 1 ) < npoints; k++ )
            {
                linear_regression( data_x + j + k * interval, data_y + j + k * interval, interval, &slice_slope, &slice_intercept );
                dx = data_x[ interval * ( k + 1 ) + j ] - data_x[ interval * k + j ];
                dy = dx * slice_slope;  
                if ( dy == dy )
                {
                    partial_sum += ( method == 0 ) ? ( ( dy > 0 ) ? dy : -dy ) : sqrt( dx * dx + dy * dy );
                    npasses++;
                }
            }
        }
        path[i] = ( partial_sum / npasses ) * ( (double)npoints / interval );
    }
    for ( i = 0; i < num_intervals; i++ )
    {
            log_intervals[i] = log( intervals[i] );
            log_path[i] = log( path[i] );
    }
    
    *correlation = linear_regression( log_intervals, log_path, num_intervals, slope, intercept );
    if ( extra_output )
    {
        printf( "path:" );
        for ( i = 0; i < num_intervals; i++ ) printf( " %g,%g", log_intervals[i], log_path[i] );
        printf( "\n" );
    }
}

#define MAX_OUTBREAKS 10
#define MAX_SAMPLES 10000

static double abs_diff( double p1, double p2 ) { if ( p1 == p2 ) return 0; return ( p1 - p2 > 0 ) ? ( p1 - p2 ) : ( p2 - p1 ); } 
static int straight_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? 1 : -1; } 
static int reversed_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? -1 : 1; } 

double averaging( double *x_coord, double *y_coord, int num_points, int method, int interval_min, int interval_max, int extra_output, int linincr, int num_parts, double *avg_slope, double *avg_intercept )
{
    int i, part_size;
    double t_slope, t_intercept, t_correlation;
    double sum_slope, sum_intercept, sum_squares;
    
    sum_slope = 0;
    sum_intercept = 0;
    sum_squares = 0;
    part_size = num_points / num_parts;
    for ( i = 0; i < num_parts; i++ )
    {
        if ( method == 0 ) 
            peng_dimension( x_coord + i * part_size, y_coord + i * part_size, part_size, 1, interval_min, interval_max, extra_output, linincr, &t_slope, &t_intercept, &t_correlation );
        else 
            higuchi_dimension( x_coord + i * part_size, y_coord + i * part_size, part_size, 1, interval_min, interval_max, extra_output, linincr, &t_slope, &t_intercept, &t_correlation );
        sum_slope += t_slope;
        sum_intercept += t_intercept;
        sum_squares += t_slope * t_slope;
    }
    *avg_slope = sum_slope / num_parts;
    *avg_intercept = sum_intercept / num_parts;
    return sqrt( sum_squares / num_parts - ( (*avg_slope) * (*avg_slope) ) );

}

/*
void logperiodic_intervals( int min_count, int num_peaks, int direction, double *peak_positions )
{
    int i;
    double peak_at_left, peak_at_right, sum_of_intervals;
    
    peak_at_left = ( direction ) ? log( min_count + 1 ) : log( min_count + num_peaks );
    peak_at_right = ( direction ) ? log( min_count + num_peaks ) : log( min_count + 1 );
    
    for ( i = 0; i < num_peaks; i++ ) 
    {
        sum_of_intervals = ( log( min_count + i + 1 ) - log( min_count + 1 ) );
        peak_positions[ i ] = ( direction ) ? 
            ( sum_of_intervals / ( peak_at_right - peak_at_left ) ) :
            ( 1 - sum_of_intervals / ( peak_at_left - peak_at_right ) );
    }
}
*/

int logperiodic_approximation( double *x_coord, double *y_coord, int num_points, int method, int interval_min, int interval_max, int extra_output, int linincr, int *first_period, int *last_period,  double *slope, double *intercept, double *correlation )
{
    int i, j, k;
    int direction;
    double lp_period, observation_period, x_pos, lp_pos, updated_x_pos;
    double *updated_x_coord;
    double t_slope, t_intercept, t_correlation, min_correlation = 0;
        
    updated_x_coord = (double*) malloc( num_points * sizeof( double ) );
    observation_period = x_coord[ num_points - 1 ] - x_coord[0];
    if ( extra_output ) printf( "Correlation\n" );
    for ( direction = 0; direction <= 1; direction++ )
    {
        if ( extra_output ) printf( "{ /*%s*/\n", ( direction == 0 ) ? "forwards" : "backwards" );
        for ( i = 2; i < 1000; i = ( i * 3 ) / 2 )
        {
            if ( extra_output ) printf( "%d: {", i );
            for ( j = 2; j < 100; j = ( j * 3 ) / 2 )
            {
                lp_period = log( i + j + 1 ) - log( i + 1 );
                for ( k = 0; k < num_points; k++ )
                {
                    x_pos = ( direction ) ? ( x_coord[k] - x_coord[0] ) : ( x_coord[ num_points - 1 ] - x_coord[k] );
                    lp_pos = log( i + 1 + j * x_pos / observation_period ) - log( i + 1 );
                    updated_x_pos = x_pos * ( lp_pos / lp_period );
                    updated_x_coord[k] = ( direction ) ? ( x_coord[0] + updated_x_pos ) : ( x_coord[ num_points - 1 ] - updated_x_pos );
                }
                
                if ( method == 0 ) 
                    peng_dimension( updated_x_coord, y_coord, num_points, 1, interval_min, interval_max, 0, linincr, &t_slope, &t_intercept, &t_correlation );
                else 
                    higuchi_dimension( updated_x_coord, y_coord, num_points, 1, interval_min, interval_max, 0, linincr, &t_slope, &t_intercept, &t_correlation );
                
                if ( extra_output ) printf( "%d:%g,", j, t_correlation );
                
                if ( t_correlation < min_correlation )
                {
                    min_correlation = t_correlation;
                    *first_period = ( direction ) ? i : i + j;
                    *last_period = ( direction ) ? i + j : i;
                    *slope = t_slope;
                    *intercept = t_intercept;
                    *correlation = t_correlation;
                    
                }
            }
            if ( extra_output ) printf( "},\n" );
        }
        if ( extra_output ) printf( "},\n" );
    }

    if ( method == 0 ) peng_dimension( x_coord, y_coord, num_points, 1, interval_min, interval_max, extra_output, linincr, &t_slope, &t_intercept, &t_correlation );
    else higuchi_dimension( x_coord, y_coord, num_points, 1, interval_min, interval_max, extra_output, linincr, &t_slope, &t_intercept, &t_correlation );

    if ( extra_output )
    {
        direction = ( *last_period > *first_period );
        i = ( direction ) ? *first_period : *last_period;
        j = ( direction ) ? *last_period - i : *first_period - i;

        lp_period = log( i + j + 1 ) - log( i + 1 );
        for ( k = 0; k < num_points; k++ )
        {
            x_pos = ( direction ) ? ( x_coord[k] - x_coord[0] ) : ( x_coord[ num_points - 1 ] - x_coord[k] );
            lp_pos = log( i + 1 + j * x_pos / observation_period ) - log( i + 1 );
            updated_x_pos = x_pos * ( lp_pos / lp_period );
            updated_x_coord[k] = ( direction ) ? ( x_coord[0] + updated_x_pos ) : ( x_coord[ num_points - 1 ] - updated_x_pos );
        }
        
        if ( method == 0 ) 
            peng_dimension( updated_x_coord, y_coord, num_points, 1, interval_min, interval_max, 1, linincr, &t_slope, &t_intercept, &t_correlation );
        else 
            higuchi_dimension( updated_x_coord, y_coord, num_points, 1, interval_min, interval_max, 1, linincr, &t_slope, &t_intercept, &t_correlation );
        
    }
    free( updated_x_coord );
    return ( t_correlation > min_correlation );
}

int option_selected( const char *key, int argc, char **argv )
{
    int i;
    for ( i = 1; i < argc; i++ ) if ( strcmp( argv[i], key ) == 0 ) return 1;
    return 0;
}

const char *option_parameter( const char *key, int argc, char **argv )
{
    int i;
    for ( i = 1; i < argc; i++ ) if ( strcmp( argv[i], key ) == 0 ) return argv[ i + 1 ];
    return 0;
}

/*
const double *array_to_compare;
int compare_func( const void *p1, const void *p2 ) { };

void sort_xy( double *x_coord, double *y_coord, int npoints )
{
    int i;
    int *order = (int*) malloc( npoints * sizeof( int ) );
    
    for ( i = 0; i < npoints; i++ ) order[i] = i;
    array_to_compare = x_coord;
}
*/

int main( int argc, char **argv )
{
    FILE *ifile = stdin;
    char buf[1024] = "";
    char *p;
    int i, maxnpoints = 0;
    int npoints = 0;
    int interval_min, interval_max, extra_output, linincr;
    double *values = 0, *xvalues = 0, *yvalues = 0;
    double *ptr;
    double slope, intercept, correlation = 0, variation;
    int first_period = 0, last_period = 0, credibility;
        
    if ( argc == 1 ) 
    { 
        printf( "Selection of algorithm:\n\n-h    \"canonical\" Higuchi f.d;\n-hm   f.d for sum of hypotenuses;\n-lh   fitting of log-periodicity\n\n-p    \"canonical\" Peng f.d\n-pm   f.d for sum of lengths\n-lp   fitting of log-periodicity\n\n-ap, -ah <n> averaging\nOptional keys:\n-xy: data in x-y representation;\n-max,-min: upper and lower limits of intervals in logarithmic dependency, in integer units;\n-v: extra output;\n-linincr: linear increment of x-coordinate in log-log dependency.\n\nstdin: input data, comma-separated.\n" ); return 1;
    }
    
    p=buf;
    while ( fread( p, 1, 1, ifile ) == 1 )
    {
        if ( *p == ',' || *p == '\n' )
        {
            if ( npoints >= maxnpoints - 1 )
            {
                maxnpoints += 1000;
                ptr = (double*) malloc( maxnpoints * sizeof( double ) );
                memcpy( ptr, values, npoints * sizeof( double ) );
                if ( values != 0 ) free( values );
                values = ptr;
            }
            *(p + 1) = 0;
            if ( strlen(buf) ) values[ npoints++ ] = atof( buf );
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
        if ( option_selected( "-xy", argc, argv ) )
        {
            npoints /= 2;
            xvalues = (double*) malloc( sizeof( double ) * npoints );
            yvalues = (double*) malloc( sizeof( double ) * npoints );
            for ( i = 0; i < npoints; i++ )
            {
                xvalues[i] = values[ 2 * i ];
                yvalues[i] = values[ 2 * i + 1 ];
            }
        }
        else
        {
            xvalues = (double*) malloc( sizeof( double ) * npoints );
            for ( i = 0; i < npoints; i++ ) xvalues[i] = i;
            yvalues = values;
        }
        extra_output = option_selected( "-v", argc, argv );
        linincr = option_selected( "-linincr", argc, argv );
        interval_min = ( option_selected( "-min", argc, argv ) ? atoi( option_parameter( "-min", argc, argv ) ) : 0 );
        interval_max = ( option_selected( "-max", argc, argv ) ? atoi( option_parameter( "-max", argc, argv ) ) : 0 );
        switch ( argv[1][1] )
        {
            case 'h':
            {
                higuchi_dimension( xvalues, yvalues, npoints, ( argv[1][2] == 'm' ) ? 1 : 0, interval_min, interval_max, extra_output, linincr, &slope, &intercept, &correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", -slope, intercept, correlation );        
                break;
            }
            case 'p':
            {
                peng_dimension( xvalues, yvalues, npoints, ( argv[1][2] == 'm' ) ? 1 : 0, interval_min, interval_max, extra_output, linincr, &slope, &intercept, &correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", -slope, intercept, correlation );        
                break;
            }
            case 'a':
            {
                variation = averaging( xvalues, yvalues, npoints, argv[1][2] == 'h', interval_min, interval_max, extra_output, linincr, atoi( argv[2] ), &slope, &intercept );
                printf( "dimension: %g\nintercept: %g\nvariation of dimension: %g\n", -slope, intercept, variation );                 
                break;
            }
            case 'l':
            {
                credibility = logperiodic_approximation( xvalues, yvalues, npoints, argv[1][2] == 'h', interval_min, interval_max, extra_output, linincr, &first_period, &last_period, &slope, &intercept, &correlation );
                printf( "credibility: %s\ndirection: %s\nperiod in beginning: %d\nperiod in end: %d\ndimension: %g\nintercept: %g\ncorrelation: %g\n", credibility ? "sufficient" : "insufficient", ( first_period < last_period ) ? "expansion from past" : "contraction to future", first_period, last_period, -slope, intercept, correlation );          
                break;
            }
            default:
                printf( "incorrect argument %s\n", argv[1] );
                break;
        }
    }
    if ( xvalues != 0 ) free( xvalues );
    if ( yvalues != 0 && yvalues != values ) free( yvalues );
    if ( values != 0 ) free( values );
    return 0;
}
