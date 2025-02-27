#define lapack_int     int
#define lapack_logical          lapack_int
typedef struct { double real, imag; } _lapack_complex_double;
#define lapack_complex_double _lapack_complex_double
#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102
#define LAPACK_WORK_MEMORY_ERROR       -1010
#define LAPACK_TRANSPOSE_MEMORY_ERROR  -1011
#define LAPACKE_malloc( size )   malloc( size )
#define LAPACKE_free( p )      free( p )
#define LAPACK_Z2INT( x ) (lapack_int)(*((double*)&x ))
#define LAPACK_DISNAN( x ) ( x != x )
#define LAPACK_ZISNAN( x ) ( LAPACK_DISNAN(*((double*)&x)) || \
                              LAPACK_DISNAN(*(((double*)&x)+1)) )

#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_

#define LAPACK_zgeev LAPACK_GLOBAL(zgeev,ZGEEV)
void LAPACK_zgeev( char* jobvl, char* jobvr, lapack_int* n,
                   lapack_complex_double* a, lapack_int* lda,
                   lapack_complex_double* w, lapack_complex_double* vl,
                   lapack_int* ldvl, lapack_complex_double* vr,
                   lapack_int* ldvr, lapack_complex_double* work,
                   lapack_int* lwork, double* rwork, lapack_int *info );

#define LAPACK_lsame LAPACK_GLOBAL(lsame,LSAME)
lapack_logical LAPACK_lsame( char* ca,  char* cb,
                              lapack_int lca, lapack_int lcb );

void LAPACKE_xerbla( const char *name, lapack_int info )
{
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        printf( "Not enough memory to allocate work array in %s\n", name );
    } else if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
        printf( "Not enough memory to transpose matrix in %s\n", name );
    } else if( info < 0 ) {
        printf( "Wrong parameter %d in %s\n", -(int) info, name );
    }
}

lapack_logical LAPACKE_lsame( char ca,  char cb )
{
    return (lapack_logical) LAPACK_lsame( &ca, &cb, 1, 1 );
}

void LAPACKE_zge_trans( int matrix_layout, lapack_int m, lapack_int n,
                        const lapack_complex_double* in, lapack_int ldin,
                        lapack_complex_double* out, lapack_int ldout )
{
    lapack_int i, j, x, y;

    if( in == NULL || out == NULL ) return;

    if( matrix_layout == LAPACK_COL_MAJOR ) {
        x = n;
        y = m;
    } else if ( matrix_layout == LAPACK_ROW_MAJOR ) {
        x = m;
        y = n;
    } else {
        /* Unknown input layout */
        return;
    }

    /* In case of incorrect m, n, ldin or ldout the function does nothing */
    for( i = 0; i < MIN( y, ldin ); i++ ) {
        for( j = 0; j < MIN( x, ldout ); j++ ) {
            out[ (size_t)i*ldout + j ] = in[ (size_t)j*ldin + i ];
        }
    }
}

lapack_logical LAPACKE_zge_nancheck( int matrix_layout, lapack_int m,
                                      lapack_int n,
                                      const lapack_complex_double *a,
                                      lapack_int lda )
{
    lapack_int i, j;

    if( a == NULL ) return (lapack_logical) 0;

    if( matrix_layout == LAPACK_COL_MAJOR ) {
        for( j = 0; j < n; j++ ) {
            for( i = 0; i < MIN( m, lda ); i++ ) {
                if( LAPACK_ZISNAN( a[i+(size_t)j*lda] ) )
                    return (lapack_logical) 1;
            }
        }
    } else if ( matrix_layout == LAPACK_ROW_MAJOR ) {
        for( i = 0; i < m; i++ ) {
            for( j = 0; j < MIN( n, lda ); j++ ) {
                if( LAPACK_ZISNAN( a[(size_t)i*lda+j] ) )
                    return (lapack_logical) 1;
            }
        }
    }
    return (lapack_logical) 0;
}

lapack_int LAPACKE_zgeev_work( int matrix_layout, char jobvl, char jobvr,
                               lapack_int n, lapack_complex_double* a,
                               lapack_int lda, lapack_complex_double* w,
                               lapack_complex_double* vl, lapack_int ldvl,
                               lapack_complex_double* vr, lapack_int ldvr,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_zgeev( &jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
                      work, &lwork, rwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,n);
        lapack_int ldvl_t = MAX(1,n);
        lapack_int ldvr_t = MAX(1,n);
        lapack_complex_double* a_t = NULL;
        lapack_complex_double* vl_t = NULL;
        lapack_complex_double* vr_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -6;
            LAPACKE_xerbla( "LAPACKE_zgeev_work", info );
            return info;
        }
        if( ldvl < n ) {
            info = -9;
            LAPACKE_xerbla( "LAPACKE_zgeev_work", info );
            return info;
        }
        if( ldvr < n ) {
            info = -11;
            LAPACKE_xerbla( "LAPACKE_zgeev_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 ) {
            LAPACK_zgeev( &jobvl, &jobvr, &n, a, &lda_t, w, vl, &ldvl_t, vr,
                          &ldvr_t, work, &lwork, rwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        if( LAPACKE_lsame( jobvl, 'v' ) ) {
            vl_t = (lapack_complex_double*)
                LAPACKE_malloc( sizeof(lapack_complex_double) *
                                ldvl_t * MAX(1,n) );
            if( vl_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_1;
            }
        }
        if( LAPACKE_lsame( jobvr, 'v' ) ) {
            vr_t = (lapack_complex_double*)
                LAPACKE_malloc( sizeof(lapack_complex_double) *
                                ldvr_t * MAX(1,n) );
            if( vr_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_2;
            }
        }
        /* Transpose input matrices */
        LAPACKE_zge_trans( matrix_layout, n, n, a, lda, a_t, lda_t );
        /* Call LAPACK function and adjust info */
        LAPACK_zgeev( &jobvl, &jobvr, &n, a_t, &lda_t, w, vl_t, &ldvl_t, vr_t,
                      &ldvr_t, work, &lwork, rwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_t, lda_t, a, lda );
        if( LAPACKE_lsame( jobvl, 'v' ) ) {
            LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, vl_t, ldvl_t, vl, ldvl );
        }
        if( LAPACKE_lsame( jobvr, 'v' ) ) {
            LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, vr_t, ldvr_t, vr, ldvr );
        }
        /* Release memory and exit */
        if( LAPACKE_lsame( jobvr, 'v' ) ) {
            LAPACKE_free( vr_t );
        }
exit_level_2:
        if( LAPACKE_lsame( jobvl, 'v' ) ) {
            LAPACKE_free( vl_t );
        }
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_zgeev_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_zgeev_work", info );
    }
    return info;
}

lapack_int LAPACKE_zgeev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* w,
                          lapack_complex_double* vl, lapack_int ldvl,
                          lapack_complex_double* vr, lapack_int ldvr )
{
    lapack_int info = 0;
    lapack_int lwork = -1;
    double* rwork = NULL;
    lapack_complex_double* work = NULL;
    lapack_complex_double work_query;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        LAPACKE_xerbla( "LAPACKE_zgeev", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    /* Optionally check input matrices for NaNs */
    if( LAPACKE_zge_nancheck( matrix_layout, n, n, a, lda ) ) {
        return -5;
    }
#endif
    /* Allocate memory for working array(s) */
    rwork = (double*)LAPACKE_malloc( sizeof(double) * MAX(1,2*n) );
    if( rwork == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    /* Query optimal working array(s) size */
    info = LAPACKE_zgeev_work( matrix_layout, jobvl, jobvr, n, a, lda, w, vl,
                               ldvl, vr, ldvr, &work_query, lwork, rwork );
    if( info != 0 ) {
        goto exit_level_1;
    }
    lwork = LAPACK_Z2INT( work_query );
    /* Allocate memory for work arrays */
    work = (lapack_complex_double*)
        LAPACKE_malloc( sizeof(lapack_complex_double) * lwork );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_1;
    }
    /* Call middle-level interface */
    info = LAPACKE_zgeev_work( matrix_layout, jobvl, jobvr, n, a, lda, w, vl,
                               ldvl, vr, ldvr, work, lwork, rwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_1:
    LAPACKE_free( rwork );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        LAPACKE_xerbla( "LAPACKE_zgeev", info );
    }
    return info;
}
