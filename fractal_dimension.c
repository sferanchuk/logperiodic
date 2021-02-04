
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

/*
void higuchi_dimension( double *data, int npoints, double *slope, double *intercept, double *correlation )
{
    int i, j, k, npasses;
    int num_logpoints, interval;
    double d, partial_sum;
    double log_path[1000], log_interval[1000];
  
    if ( npoints > 100 ) { num_logpoints = (int) ( log( npoints ) / log( 2 ) ); for ( i = 1; i < num_logpoints; i++ ) log_interval[i] = i * log( 2 ); }
    else { num_logpoints = npoints - 3; for ( i = 0; i < num_logpoints; i++ ) log_interval[i] = log( i + 2 ); }
    
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( log_interval[i] );
        partial_sum = 0;
        npasses = 0;
        for ( j = 0; j < interval; j++ )
        {
            for ( k = 0; j + interval * ( k + 1 ) < npoints; k++ )
            {
                d = data[ interval * ( k + 1 ) + j ] - data[ interval * k + j ];
                partial_sum += ( d > 0) ? d : -d;
                npasses++;
            }
        }
        log_path[i] = -log( ( partial_sum / npasses ) * ( (double)npoints / interval ) );
    }
    *correlation = linear_regression( log_interval, log_path, num_logpoints, slope, intercept );
 }
*/

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
*/

int extra_output = 0;

void higuchi_dimension( double *data_x, double *data_y, int npoints, int method, double *slope, double *intercept, double *correlation )
{
    int i, j, k, npasses;
    int num_logpoints, interval;
    double dx, dy, partial_sum;
    double log_path[1000], log_interval[1000];
  
    if ( npoints > 100 ) 
    { 
        num_logpoints = (int) ( log( npoints ) / log( 2 ) ); 
        for ( i = 0; i < num_logpoints; i++ ) log_interval[i] = i * log( 2 ); 
    }
    else 
    { 
        num_logpoints = npoints - 3; 
        for ( i = 0; i < num_logpoints; i++ ) log_interval[i] = log( i + 2 ); 
    }
    
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( log_interval[i] );
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
        log_path[i] = -log( ( partial_sum / npasses ) * ( (double)npoints / interval ) );
    }
    *correlation = linear_regression( log_interval, log_path, num_logpoints, slope, intercept );
    if ( extra_output )
    {
        printf( "path:" );
        for ( i = 0; i < num_logpoints; i++ ) printf( " %g,%g", log_interval[i], log_path[i] );
        printf( "\n" );
    }
}


void peng_dimension( double *data_x, double *data_y, int npoints, double *slope, double *intercept, double *correlation )
{
    int i, j, k, npasses;
    int num_logpoints, interval;
    double slice_slope, slice_intercept, dx, dy, partial_sum;
    double log_path[1000], log_interval[1000];
  
    if ( npoints > 100 ) 
    { 
        num_logpoints = (int) ( log( npoints ) / log( 2 ) ); 
        for ( i = 0; i < num_logpoints; i++ ) log_interval[i] = ( i + 1 ) * log( 2 ); 
    }
    else 
    { 
        num_logpoints = npoints - 3; 
        for ( i = 0; i < num_logpoints; i++ ) log_interval[i] = log( i + 3 ); 
    }
    
    for ( i = 0; i < num_logpoints; i++ )
    {
        interval = (int) exp( log_interval[i] );
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
                    partial_sum += ( dy > 0) ? dy : -dy;
                npasses++;
            }
        }
        log_path[i] = -log( ( partial_sum / npasses ) * ( (double)npoints / interval ) );
    }
    *correlation = linear_regression( log_interval, log_path, num_logpoints, slope, intercept );
    if ( extra_output )
    {
        printf( "path:" );
        for ( i = 0; i < num_logpoints; i++ ) printf( " %g,%g", log_interval[i], log_path[i] );
        printf( "\n" );
    }
}

#define MAX_OUTBREAKS 10
#define MAX_SAMPLES 10000

static double abs_diff( double p1, double p2 ) { if ( p1 == p2 ) return 0; return ( p1 - p2 > 0 ) ? ( p1 - p2 ) : ( p2 - p1 ); } 
static int straight_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? 1 : -1; } 
static int reversed_compare( const void *p1, const void *p2 ) { if ( *(double*)p1 == *(double*)p2 ) return 0; return ( *(double*)p1 - *(double*)p2 > 0 ) ? -1 : 1; } 

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

int logperiodic_approximation( double *x_coord, double *y_coord, int num_points, int method, int *first_period, int *last_period, double *correlation )
{
    int i, j, k;
    int direction;
    double peak_at_left, peak_at_right, sum_of_intervals;
    double *updated_x_coord;
    double t_slope, t_intercept, t_correlation, max_correlation;
        
    updated_x_coord = (double*) malloc( num_points * sizeof( double ) );
 
    max_correlation = 0;
    
    for ( i = 2; i < 1000; i = ( i * 3 ) / 2 )
    for ( j = 2; j < 100; j = ( j * 3 ) / 2 )
    for ( direction = 0; direction <= 1; direction++ )
    {
        
        peak_at_left = ( direction ) ? log( i + 1 ) : log( i + j );
        peak_at_right = ( direction ) ? log( i + j ) : log( i + 1 );
        for ( k = 0; k < i; k++ )
        {
                sum_of_intervals =  log( i + (double)( k * j ) / num_points ) - log( i + 1 );
                updated_x_coord[k] = ( k ) ?
                    ( x_coord[0] + ( x_coord[k] - x_coord[0] ) * sum_of_intervals / ( peak_at_right - peak_at_left ) ) :
                    ( x_coord[ num_points - 1 ] - ( x_coord[ num_points - 1 ] - x_coord[k] ) * sum_of_intervals / ( peak_at_left - peak_at_right ) );
        }
        if ( method == 0 ) peng_dimension( updated_x_coord, y_coord, num_points, &t_slope, &t_intercept, &t_correlation );
        else higuchi_dimension( updated_x_coord, y_coord, num_points, 1, &t_slope, &t_intercept, &t_correlation );
        /*if ( extra_output ) printf( "%d %d %g %g\n", i, j, t_slope, t_correlation );*/
        if ( t_correlation > max_correlation )
        {
            max_correlation = t_correlation;
            *first_period = ( direction ) ? i : i + j;
            *last_period = ( direction ) ? i + j : i;
            *correlation = max_correlation;
            
        }
    }
    free( updated_x_coord );

    if ( method == 0 ) peng_dimension( x_coord, y_coord, num_points, &t_slope, &t_intercept, &t_correlation );
    else higuchi_dimension( x_coord, y_coord, num_points, 1, &t_slope, &t_intercept, &t_correlation );
    return ( t_correlation < max_correlation );
}

    /*
    int i, j, k, num_samples, num_outbreaks, max_num_periods;
    double sum, sum_sq, d, vmin, slope, intercept;
    double eps, global_mean, local_mean, global_sd, local_sd, outbreak_bound;
    double log_distr[ MAX_SAMPLES ], log_count[ MAX_SAMPLES ];
    int indices[ MAX_OUTBREAKS] ;
    double outbreak_x[ MAX_OUTBREAKS ], outbreak_y[ MAX_OUTBREAKS ], *ptr_x, *ptr_y;
    double peak_positions[ MAX_OUTBREAKS ];

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
    max_num_periods = 1000;


    vmin = -1;
    for ( i = 1; i < max_num_periods; i *= 2 )
    {
        for ( j = 0; j < 2; j++ )
        {
            logperiodic_intervals( i, num_outbreaks, j, peak_positions );

            d = linear_regression( peak_positions, outbreak_x, num_outbreaks, &slope, &intercept );
 
            if ( vmin < 0 || vmin > d )
            {
                vmin = d;
                *direction = j;
                *first_period = i;   
                *rate = slope;
                *phase = intercept;
                *correlation = d;
            }
        }
    }
    
 
    if ( print_datapoints ) 
    {
        printf( "outbreaks_x: " );
        for ( j = 0; j < num_outbreaks; j++ ) 
            printf( "%7g ", outbreak_x[j] );
        printf( "\n" );
        printf( "outbreaks_y: " );
        for ( i = 0; i < num_outbreaks; i++ ) 
        for ( j = 0; j < num_outbreaks; j++ )
            if ( data_x[ indices[j] ] == outbreak_x[i] ) printf( "%7g ", data_y[ indices[j] ] ); 
        printf( "\n" );
        printf( "approximated periods: " );
        logperiodic_intervals( *first_period, num_outbreaks, *direction, peak_positions );
        for ( j = 0; j < num_outbreaks; j++ ) 
            printf( "%7g ", peak_positions[ j ] * (*rate) + (*phase) );
        printf( "\n" );
    }
}
*/

    
int option_selected( const char *key, int argc, char **argv )
{
    int i;
    for ( i = 1; i < argc; i++ ) if ( strcmp( argv[i], key ) == 0 ) return 1;
    return 0;
}

int main( int argc, char **argv )
{
    FILE *ifile = stdin;
    char buf[1024] = "";
    char *p;
    int i, maxnpoints = 0;
    int npoints = 0;
    double *values = 0, *xvalues = 0, *yvalues = 0;
    double *ptr;
    double correlation, dimension, intercept;
    int first_period, last_period, credibility;
        
    if ( argc == 1 ) 
    { 
        printf( "-h: \"canonical\" higuchi f.d; -hm: higuchi-style f.d. with variable x coord;\n-p: peng f.d; -c,-d: log-periodic approximation;\noptional keys: -xy : data in x-y representation; -v extra output.\nstdin: input data, comma-separated, in one line.\n" ); return 1;
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
        switch ( argv[1][1] )
        {
            case 'h':
            {
                higuchi_dimension( xvalues, yvalues, npoints, ( argv[1][2] == 'm' ) ? 1 : 0, &dimension, &intercept, &correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", dimension, intercept, correlation );        
                break;
            }
            case 'p':
            {
                peng_dimension( xvalues, yvalues, npoints,  &dimension, &intercept, &correlation );
                printf( "dimension: %g\nintercept: %g\ncorrelation: %g\n", dimension, intercept, correlation );        
                break;
            }
            case 'c':
            case 'd':
            {
                credibility = logperiodic_approximation( xvalues, yvalues, npoints, argv[1][1] == 'c', &first_period, &last_period, &correlation );
                printf( "credibility: %s\ndirection: %s\nperiod in beginning: %d\nperiod in end: %d\ncorrelation: %g\n", credibility ? "sufficient" : "insufficient", ( first_period < last_period ) ? "expansion from past" : "contraction to future", first_period, last_period, correlation );          
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
