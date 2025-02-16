!---------------------------------------------------------------
! ocean2cmaq:
!       Prepare OCEAN_1 file for CMAQ run.
!
! author: 
!       Ramiro A. Espada. 
!       Lakes Environmental Software
!       February 2025.
!---------------------------------------------------------------
program ocean2cmaq

  use netcdf

  implicit none
  !Parameters:
  real, parameter    :: R_EARTH = 6370000.
  real, parameter    :: PI = 3.141592653589793
  real, parameter    :: RAD2DEG = 180./PI, DEG2RAD = PI/180.

  !Objects/Strucs:
  type proj_type
     character(16)    :: pName                                    !Projection name
     integer          :: typ                                      !Integer code for projection TYPE (2=lcc, 6=stere, 7=merc)
     real             :: alp,bet,gam,xcent,ycent                  !proj parameters.
     real             :: p1,p2,p3,p4                              !extra parameters to speed up calculation once p%typ is defined.
  end type proj_type

  type grid_type
      character(12)   :: gName                                    !grid-name
      integer         :: nx,ny,nz                                 !number of cells in x-y direction (ncols, nrows, nlevs)
      real            :: dx,dy                                    !x-y cell dimension (x_cell, y_cell)
      real            :: xmin,ymin,xmax,ymax,xc,yc
      real            :: lonmin,latmin,lonmax,latmax
  end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  character(4)   :: YYYY
  character(2)   :: MM,DD
  character(254) :: outFile

  integer :: status,iostat,ncid,varid,dimid
  integer :: i,j,k
  !coordinates
  !real, allocatable, save  :: longitude(:,:), latitude(:,:)
  real, allocatable, save  :: opn(:,:), surf(:,:), dms(:,:), chlo(:,:)
  character(len=16) :: var_list(4),var_unit(4),var_dtyp(4)

  character(len=17) :: start_date!,end_date
  character(200)    :: griddesc_file,gridname,open_file,surf_file,chlo_file,dms_file
  namelist/control/start_date,griddesc_file,gridname,open_file,surf_file,chlo_file,dms_file

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'ocean2cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file,gridname, proj, grid)
   
  print*,grid%gName
  print*,grid%dx, grid%dy, grid%nx,grid%ny
  print*, grid%xmin,grid%ymin,grid%xmax,grid%ymax,grid%xc,grid%yc  
  print*, grid%lonmin,grid%latmin,grid%lonmax,grid%latmax

  allocate( opn(grid%nx, grid%ny))
  allocate(surf(grid%nx, grid%ny))
  allocate(chlo(grid%nx, grid%ny))
  allocate( dms(grid%nx, grid%ny))

  YYYY=start_date(1:4)
  MM=start_date(6:8)
  DD=start_date(10:12)
  outFile="ocean_"//MM//"_"//trim(grid%gName)//".nc"
  !(1) Create NetCDF "ocean_MM_GRIDNAME.nc" with IOAPI standars:

  var_list(1)="OPEN"; var_unit(1)="      " ; var_dtyp(1)="FLOAT" 
  var_list(2)="SURF"; var_unit(2)="      " ; var_dtyp(2)="FLOAT"
  var_list(3)="CHLO"; var_unit(3)="nM    " ; var_dtyp(3)="FLOAT"
  var_list(4)="DMS" ; var_unit(4)="mg m-3" ; var_dtyp(4)="FLOAT"
  call createNetCDF(trim(outFile),proj,grid, var_list, var_unit, var_list, var_dtyp)


  print*, "copy and paste"
  !(2) Put everything insde the ocean_<MM>_<GRIDNAME>.nc file.
  !(2a) Read open and copy it inside ocean_MM_GRIDNAME.nc
  print*, "OPEN..",open_file
  call check(nf90_open(trim(open_file), nf90_nowrite, ncid ))
    call check( nf90_inq_varid(ncid,"OPEN", varid        ))
    call check( nf90_get_var(ncid, varid , opn           ))
  call check(nf90_close(ncid))

  !(2b) Read SURF and copy it inside ocean_MM_GRIDNAME.nc
  print*, "SURF..",surf_file
  call check(nf90_open(trim(surf_file), nf90_nowrite, ncid ))
    call check( nf90_inq_varid(ncid,"SURF", varid        ))
    call check( nf90_get_var(ncid, varid , surf          ))
  call check(nf90_close(ncid))
  !(2c) Read CHLO and copy it inside ocean_MM_GRIDNAME.nc

  print*, "CHLO..",chlo_file
  call check(nf90_open(trim(chlo_file), nf90_nowrite, ncid ))
    call check( nf90_inq_varid(ncid,"CHLO", varid        ))
    call check( nf90_get_var(ncid, varid , chlo          ))
  call check(nf90_close(ncid))
  !(2d) Read DMS  and copy it inside ocean_MM_GRIDNAME.nc

  print*, "DMS..",dms_file
  call check(nf90_open(trim( dms_file), nf90_nowrite, ncid ))
    call check( nf90_inq_varid(ncid,"DMS" , varid        ))
    call check( nf90_get_var(ncid, varid , dms           ))
  call check(nf90_close(ncid))


  !WRITE EVERYTHING ON OUTPUT FILE:
  print*,"write", outFile
  call check(nf90_open(trim(outFile), nf90_write, ncid ))
       call check(nf90_inq_varid(ncid,"OPEN"    ,varid)); call check(nf90_put_var(ncid, varid,  OPN(:,:) ))
       call check(nf90_inq_varid(ncid,"SURF"    ,varid)); call check(nf90_put_var(ncid, varid, SURF(:,:) ))
       call check(nf90_inq_varid(ncid,"CHLO"    ,varid)); call check(nf90_put_var(ncid, varid, CHLO(:,:) ))
       call check(nf90_inq_varid(ncid,"DMS"     ,varid)); call check(nf90_put_var(ncid, varid,  DMS(:,:) ))

       call check(nf90_inq_varid(ncid, "TFLAG"  ,varid));
       call check(nf90_put_var(ncid, varid, reshape([0, 0, 0, 0, 0, 0, 0, 0], [2,4,1])))
  call check(nf90_close(ncid))

print*, "========================================="
print*, " ocean2cmaq: Completed successfully"
print*, "========================================="

contains
!UTILS *****************************************************************
 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 
subroutine createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)
    implicit none
    type(grid_type) , intent(in) :: g
    type(proj_type) , intent(in) :: p
    character(len=254), intent(in) :: outFile
    character(len=16) :: var_list(:),var_unit(:),var_dtyp(:)
    character(len=25) :: var_desc(:)
    integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
    integer :: k
    integer :: nvars
    character(6000) :: var_list_string
 
    nvars=size(var_list)
    write(var_list_string,*) var_list                !este es un global attr importante.
    
    call check(nf90_create(outFile, NF90_CLOBBER, ncid))
        ! Defino dimensiones
        call check(nf90_def_dim(ncid, "TSTEP"    , 1      , tstep_dim_id     ))
        call check(nf90_def_dim(ncid, "DATE-TIME", 2      , date_time_dim_id ))
        call check(nf90_def_dim(ncid, "COL"      , g%nx   , col_dim_id       ))
        call check(nf90_def_dim(ncid, "ROW"      , g%ny   , row_dim_id       ))
        call check(nf90_def_dim(ncid, "LAY"      , 1      , lay_dim_id       ))
        call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id       ))
        !Defino variables
        call check(nf90_def_var(ncid,"TFLAG",NF90_INT      , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
        call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
        call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))

        do k=1, nvars
          if ( trim(var_dtyp(k)) == "INT" ) then
                call check(nf90_def_var(ncid, trim(var_list(k)) , NF90_INT  , [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) 
          else if ( trim(var_dtyp(k)) == "FLOAT" ) then
                call check(nf90_def_var(ncid, trim(var_list(k)) , NF90_FLOAT, [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) 
          endif
          call check(nf90_put_att(ncid, var_id,"long_name",      var_list(k)  ))
          call check(nf90_put_att(ncid, var_id,"units"    , trim(var_unit(k)) ))
          call check(nf90_put_att(ncid, var_id,"var_desc" , trim(var_desc(k)) ))
        end do
        ! Defino attributos
        call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
        call check(nf90_put_att(ncid, nf90_global,"EXEC_ID", "????????????????"   ))
        call check(nf90_put_att(ncid, nf90_global,"FTYPE"  , 1                    ))
        call check(nf90_put_att(ncid, nf90_global,"SDATE"  , 0000000              ))!stat_date (int)
        call check(nf90_put_att(ncid, nf90_global,"STIME"  , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"WDATE"  , 0000000              ))
        call check(nf90_put_att(ncid, nf90_global,"WTIME"  , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"CDATE"  , 0000000              ))
        call check(nf90_put_att(ncid, nf90_global,"CTIME"  , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"TSTEP"  , 10000                ))
        call check(nf90_put_att(ncid, nf90_global,"NTHIK"  , 1                    ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"NCOLS"  , g%nx                 ))
        call check(nf90_put_att(ncid, nf90_global,"NROWS"  , g%ny                 ))
        call check(nf90_put_att(ncid, nf90_global,"NLAYS"  , 1                    ))!grid%nz
        call check(nf90_put_att(ncid, nf90_global,"NVARS"  , nvars                ))
        call check(nf90_put_att(ncid, nf90_global,"GDTYP"  , p%typ                ))
        call check(nf90_put_att(ncid, nf90_global,"P_ALP"  , p%alp                ))
        call check(nf90_put_att(ncid, nf90_global,"P_BET"  , p%bet                ))
        call check(nf90_put_att(ncid, nf90_global,"P_GAM"  , p%gam                ))
        call check(nf90_put_att(ncid, nf90_global,"XCENT"  , p%xcent              ))
        call check(nf90_put_att(ncid, nf90_global,"YCENT"  , p%ycent              ))
        call check(nf90_put_att(ncid, nf90_global,"XORIG"  , g%xmin               ))
        call check(nf90_put_att(ncid, nf90_global,"YORIG"  , g%ymin               ))
        call check(nf90_put_att(ncid, nf90_global,"XCELL"  , g%dx                 ))
        call check(nf90_put_att(ncid, nf90_global,"YCELL"  , g%dy                 ))
        call check(nf90_put_att(ncid, nf90_global,"VGTYP"  , -9999                ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"VGTOP"  , 0.                   ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"VGLVLS" , [0., 0.]             ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"GDNAM"  , g%gName              ))
        call check(nf90_put_att(ncid, nf90_global,"UPNAM"  , "ocean2cmaq.exe"     ))!no sé que es.
        call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string)))
        call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "OCEAN_1 input file"   ))
        call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
     call check(nf90_enddef(ncid))
     !End NetCDF define mode
 end subroutine createNetCDF
!***********************************************************************
!GRIDDESC **************************************************************
 subroutine read_GRIDDESC(griddescFile,gridName, p, g)
  implicit none
  character(200),intent(in) :: griddescFile
  character(*) ,intent(in)  :: gridName
  type(proj_type), intent(inout) :: p
  type(grid_type), intent(inout) :: g
  character(20) :: row
  integer :: iostat
  iostat=0
  open(unit=2,file=griddescFile,status='old',action='read',access='sequential')
  do while(iostat == 0)  !loop por cada fila
     read(2,*,iostat=iostat) row
     if ( trim(row) == trim(gridname)) then
       g%gName=row
       read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny !projName xorig yorig xcell ycell nrows ncols
       rewind(2)
     endif
     if (trim(row) == trim(p%pName)) then
       read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent  
       iostat=1
     endif
  enddo
  close(2)

  call set_additional_proj_params(p)
  call set_additional_grid_params(p,g)
 end subroutine

 subroutine set_additional_proj_params(p)
    implicit none
    type(proj_type) ,intent(inout) :: p

    if ( p%typ == 1 ) then    !latlon                          
       print*, "Warning: Latlon coordinate system! (Not tested)."

    else if ( p%typ == 2 ) then       !lambert conformal conic:        
       if ( ABS(p%alp - p%bet) > 0.1 ) then  !secant proj case
          p%p2=     LOG( COS(p%alp           *deg2rad )/ COS(p%bet           *deg2rad)   )
          p%p2=p%p2/LOG( TAN((45.0+0.5*p%bet)*deg2rad )/ TAN((45.0+0.5*p%alp)*deg2rad)   ) !n
        else                                 !tangent proj case
          p%p2=SIN(p%alp*deg2rad) !n
       endif
       p%p3=R_EARTH*(COS(p%alp*deg2rad)*TAN((45+0.5*p%alp)*deg2rad)**p%p2)*(1/p%p2)  !F
       p%p1=p%p3/(TAN((45 + 0.5*p%ycent)*deg2rad)**p%p2)                             !rho0 

    else if ( p%typ == 6 ) then  !polar secant stereographic
       print*, "Todavia no desarrollado soporte para proyeccion polar stereografica (ptype=",p%typ,")."; stop

    else if ( p%typ == 7 ) then  !equatorial mercator
       p%p1=COS(p%alp*deg2rad)   !k0

    else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
    end if
 end subroutine

 subroutine set_additional_grid_params(p,g)
    implicit none
    type(proj_type) ,intent(inout) :: p
    type(grid_type) ,intent(inout) :: g
    real :: latmin,lonmin,latmax,lonmax
    !Obtener coordenadas del centro de la grilla, min y max:
    g%xc=0.0;g%yc=0.0;g%xmax=g%xmin+g%dx*g%nx; g%ymax=g%ymin+g%dy*g%ny

    !calculo minimos y maximos de latlon 
    !   (ojo! Dado que son transf no-lineales no corresponden necesariamente a los vertices)
    call xy2ll(p,g%xmin,g%ymin,g%lonmin,g%latmin)       !lower-left
    call xy2ll(p,g%xmax,g%ymax,g%lonmax,g%latmax)       !upper-right
 
    !latmin
    call xy2ll(p,g%xmin+g%dx*g%nx*0.5, g%ymin,lonmin,latmin)
    g%latmin=min(g%latmin,latmin)
    !latmax
    call xy2ll(p,g%xmax-g%dx*g%nx*0.5,g%ymax ,lonmax,latmax)
    g%latmax=max(g%latmax,latmax)

    !lonmin   
    call xy2ll(p,g%xmin,g%ymin+g%dy*g%ny*0.5,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)
    !         
    call xy2ll(p,g%xmin,g%ymax              ,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)

    !lonmax
    call xy2ll(p,g%xmax,g%ymax-g%dy*g%ny*0.5,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)
    !      
    call xy2ll(p,g%xmax,g%ymin              ,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)

 end subroutine
!***********************************************************************
!PROJ     **************************************************************
!COORDINATE TRANSFORMATION FUNCTIONS:======================================
 subroutine xy2ll(p,x,y,lon,lat)
     implicit none                            
     type(proj_type) ,intent(in) :: p
     real, intent(in)   :: x,y
     real, intent(inout):: lon,lat
 
     if      ( p%typ == 1 ) then  !lat lon                 
        lon=x;  lat=y
     else if ( p%typ == 2 ) then  !Lambert Conformal Conic:
        call xy2ll_lcc(p,x,y,lon,lat)
     else if ( p%typ == 6 ) then  !polar secant stereographic
        call xy2ll_stere(p,x,y,lat,lon)
     else if ( p%typ == 7 ) then  !equatorial mercator
        call xy2ll_merc(p,x,y,lon,lat)
     else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
     end if
 end subroutine
 subroutine ll2xy(p,lon,lat,x,y)
       implicit none                            
       type(proj_type) ,intent(in) :: p
       real, intent(in):: lon,lat
       real, intent(inout)   :: x,y
 
       if ( p%typ == 1 ) then        !latlon
           x=lon;y=lat               !no transformation needed
       else if ( p%typ == 2 ) then  !Lambert Conformal Conic:
          call ll2xy_lcc(p,lon,lat,x,y)
       else if ( p%typ == 6 ) then  !Polar Secant Stereographic
          call ll2xy_stere(p,lon,lat,x,y)
       else if ( p%typ == 7 ) then  !Equatorial Mercator
          call ll2xy_merc(p,lon,lat,x,y)
       else
          print*, "codigo de proyección invalido:",p%typ,"."; stop
       end if                             
 end subroutine

 !--------------------------------------------------------------------------
 !LAMBERT CONFORMAL CONIC:
 subroutine xy2ll_lcc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: n,F,rho0,rho,theta
   
   rho0=p%p1
   n=p%p2
   F=p%p3
   
   theta=ATAN(x/(rho0-y))*rad2deg
   rho=SIGN(1.0,n) * SQRT( x*x + (rho0-y)*(rho0-y))
   
   lon=p%gam+theta/n
   lat=2.0 * ATAN( (F/rho)**(1/n) )*rad2deg - 90.0 
 end subroutine
 subroutine ll2xy_lcc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: n,F,rho0,rho,dlon

   !interm params:
   rho0=p%p1
   n=p%p2
   F=p%p3

   rho=F/(TAN((45.0 + 0.5*lat)*deg2rad)**n)
   dlon=lon-p%gam
   !
   x=     rho*SIN(n*dlon*deg2rad )
   y=rho0-rho*COS(n*dlon*deg2rad )
 end subroutine
 !--------------------------------------------------------------------------
 !MERCATOR                
 subroutine xy2ll_merc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k0R,phi
 
   k0R=p%p1*R_EARTH     !es una longitud (k0*R_EARTH)
   phi=y/k0R            !es un angulo
   
   lon=p%gam + x/k0R*rad2deg
   lat=90.0-2*ATAN( EXP(-phi) )*rad2deg
 end subroutine
 subroutine ll2xy_merc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k0,lon0

   k0=p%p1              !adminesional
   lon0=p%gam           !es un angulo

   x=k0*R_EARTH*(lon-lon0)*deg2rad
   y=k0*R_EARTH*LOG(TAN((45.0+0.5*lat)*deg2rad))
 end subroutine
!--------------------------------------------------------------------------
 !POLAR STEREOGRAPHIC     
 subroutine xy2ll_stere(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k,rho

   stop 'Stereographic proyection not yet tested!'
   rho = sqrt(x*x+y*y)
   k = 2.0*ATAN( rho/2.0/R_EARTH )

   lat =         ASIN(   COS(k)*SIN(p%alp*deg2rad) + y*SIN(k)*COS(p%alp*deg2rad)/rho )               * rad2deg
   lon = p%gam + ATAN( x*SIN(k)  / ( rho*COS(p%alp*deg2rad)*COS(k) - y*SIN(p%alp*deg2rad)*SIN(k) ) ) * rad2deg

 end subroutine
 subroutine ll2xy_stere(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k!,hemi

   stop 'Stereographic proyection not yet tested!'
   !hemi=SIGN(1.0,p%alp)
   k = 2.0*R_EARTH / (1 + SIN(p%alp*deg2rad)*SIN(lat*deg2rad) + COS(p%alp*deg2rad)*COS(lat*deg2rad)*COS( (lon-p%gam)*deg2rad ))

   x = k *   COS( lat *deg2rad) * SIN( (lon - p%gam)*deg2rad )
   y = k * ( COS(p%alp*deg2rad) * SIN(  lat         *deg2rad ) - SIN(p%alp*deg2rad)*COS(lat*deg2rad)*COS((lon-p%gam)*deg2rad) )
 end subroutine
!!END COORDINATE TRANFORMATION FUNCTIONS====================================
end program
