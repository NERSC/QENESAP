  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from PH/ph.f90  
  !-----------------------------------------------------------------------
  PROGRAM epw
  !! author: Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !! version: v4.0
  !! license: GNU
  !! summary: EPW main driver 
  !!  
  !! This is the main EPW driver which sets the phases on the wavefunctions,
  !! calls [[wann_run]] and [[elphon_shuffle_wrap]]
  !!
  !! @Note
  !! 8/14/08 lnscf is unnecessary, as is nqs,iq_start
  !!
  USE io_global,       ONLY : stdout
#ifdef __PARA
  USE mp,              ONLY : mp_bcast, mp_barrier
  USE mp_world,        ONLY : mpime  
  USE mp_global,       ONLY : my_pool_id, mp_startup, ionode_id, mp_global_end
#endif
  USE control_flags,   ONLY : gamma_only
  USE control_epw,     ONLY : wannierize
  USE global_version,  ONLY : version_number
  USE epwcom,          ONLY : filukk, eliashberg, ep_coupling
  USE environment,     ONLY : environment_start
  USE elph2,           ONLY : elph 
  ! Flag to perform an electron-phonon calculation. If .true. 
  ! the code will enter in [[elphon_shuffle_wrap]]
  !
  !
  implicit none
  !
  CHARACTER (LEN=12)   :: code = 'EPW'
  !! Name of the program
  !
  version_number = '4.0.0'
  !
  !
  CALL init_clocks( .TRUE. )
  !
  CALL start_clock( 'EPW' )
  ! 
  !
  gamma_only = .FALSE.
  !
#ifdef __PARA
  CALL mp_startup()
#endif
  !
  ! Display the logo
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
write(stdout,'(a)') "                                                                                      "
write(stdout,'(a)') "                                       ``:oss/                                        "
write(stdout,'(a)') "                           `.+s+.     .+ys--yh+     `./ss+.                           "
write(stdout,'(a)') "                          -sh//yy+`   +yy   +yy    -+h+-oyy                           "
write(stdout,'(a)') "                          -yh- .oyy/.-sh.   .syo-.:sy-  /yh                           "
write(stdout,'(a)') "                 `.-.`    `yh+   -oyyyo.     `/syys:    oys      `.`                  "
write(stdout,'(a)') "               `/+ssys+-` `sh+      `                   oys`   .:osyo`                "
write(stdout,'(a)') "               -yh- ./syyooyo`                          .sys+/oyo--yh/                "
write(stdout,'(a)') "               `yy+    .-:-.                             `-/+/:`  -sh-                "
write(stdout,'(a)') "                /yh.                                              oys                 "
write(stdout,'(a)') "          ``..---hho---------`   .---------..`      `.-----.`    -hd+---.             "
write(stdout,'(a)') "       `./osmNMMMMMMMMMMMMMMMs. +NNMMMMMMMMNNmh+.   yNMMMMMNm-  oNMMMMMNmo++:`        "
write(stdout,'(a)') "       +sy--/sdMMMhyyyyyyyNMMh- .oyNMMmyyyyyhNMMm+` -yMMMdyyo:` .oyyNMMNhs+syy`       "
write(stdout,'(a)') "       -yy/   /MMM+.`-+/``mMMy-   `mMMh:`````.dMMN:` `MMMy-`-dhhy```mMMy:``+hs        "
write(stdout,'(a)') "        -yy+` /MMMo:-mMM+`-oo/.    mMMh:     `dMMN/`  dMMm:`dMMMMy..MMMo-.+yo`        "
write(stdout,'(a)') "         .sys`/MMMMNNMMMs-         mMMmyooooymMMNo:   oMMM/sMMMMMM++MMN//oh:          "
write(stdout,'(a)') "          `sh+/MMMhyyMMMs- `-`     mMMMMMMMMMNmy+-`   -MMMhMMMsmMMmdMMd/yy+           "
write(stdout,'(a)') "    `-/+++oyy-/MMM+.`/hh/.`mNm:`   mMMd+/////:-.`      NMMMMMd/:NMMMMMy:/yyo/:.`      "
write(stdout,'(a)') "   +os+//:-..-oMMMo:--:::-/MMMo. .-mMMd+---`           hMMMMN+. oMMMMMo. `-+osyso:`   "
write(stdout,'(a)') "   syo     `mNMMMMMNNNNNNNNMMMo.oNNMMMMMNNNN:`         +MMMMs:`  dMMMN/`     ``:syo   "
write(stdout,'(a)') "   /yh`     :syyyyyyyyyyyyyyyy+.`+syyyyyyyyo:`         .oyys:`   .oyys:`        +yh   "
write(stdout,'(a)') "   -yh-        ````````````````    `````````              ``        ``          oys   "
write(stdout,'(a)') "   -+h/------------------------::::::::://////++++++++++++++++++++++///////::::/yd:   "
write(stdout,'(a)') "   shdddddddddddddddddddddddddddddhhhhhhhhyyyyyssssssssssssssssyyyyyyyhhhhhhhddddh`   "
write(stdout,'(a)') "                                                                                      "
write(stdout,'(a)') "  S. Ponce, E. R. Margine, C. Verdi, and F. Giustino, http://arxiv.org/abs/1604.03525 "
write(stdout,'(a)') "                                                                                      "
#ifdef __PARA
  ENDIF
#endif  
  !
  CALL environment_start ( code )
  !
  IF ( ep_coupling ) & 
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  !
  ! Read in the input file
  !
  CALL epw_readin
  !
  CALL allocate_epwq
  CALL epw_setup
  !
  !  Print run info to stdout
  !
  CALL epw_summary
  !
  IF ( ep_coupling ) THEN 
     !
     CALL openfilepw
     !
     CALL print_clock( 'EPW' )
     !
     CALL epw_init(.true.)
     !
     CALL print_clock( 'EPW' )
     !
     !  Generates the perturbation matrix which fixes the gauge of 
     !  the calculated wavefunctions
     !
     CALL setphases_wrap
     !
     IF (wannierize) THEN
        !
        !  Create U(k, k') localization matrix 
        !      
        CALL wann_run
     ELSE
        !
        ! Read Wannier matrix from a previous run
        !
        WRITE(stdout,'(/,5x,a,/,3a,/,5x,a,/)') repeat('-',67), '     Using ', &
             trim(filukk) , ' from disk', repeat('-',67) 
     ENDIF
     !
     IF ( elph ) THEN
        !
        CALL dvanqq2()
        !
        CALL elphon_shuffle_wrap()
        !
    ENDIF
    !
    ! ... cleanup of the variables
    !
    CALL clean_pw( .FALSE. )
    CALL deallocate_epw
    !
    ! ... Close the files
    !
    CALL close_epw()
    !
  ENDIF
  !
  IF ( eliashberg ) THEN
     CALL eliashberg_eqs()
  ENDIF
  !
  ! ... Print statistics and exit gracefully    
  !
  CALL stop_epw( .TRUE. )
  !
  STOP
  !
  END PROGRAM epw
