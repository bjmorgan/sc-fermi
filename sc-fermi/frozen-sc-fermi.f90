program frozen_fermi
!
! Program to determine the self consistent Fermi level in a 
! semiconductor (or metal) by considering the formation energy
! of a set of defects in different charge states and the DOS
! of the perfect system.
!
! The relevant equations are Eq. (2-4) of Yang et al., Scientific
! Reports, 5:16977 (2015).
!
!*************************
! UPDATE - equations in that reference are wrong - see labbook
! May 2017
!*************************
!
! This version is identical to sc-fermi, but reads in a defect
! concentration as 'frozen in', then computes the self-consistent
! Fermi level given this frozen concentration.
!
! Given the list of 'original defects', one can enter a total 
! fixed concentration for a particular defect, and one can 
! fix particular charge states. In addition, one can add extra
! defects that weren't part of the list of 'original defects' but
! they must be given with a fixed concentration for a particular 
! charge state.
!
! We define the Gibbs factors as B_q = g_q exp(-E_f^q/kT), 
! where g_q is the degeneracy factor, E_f^q the formation energy
! and kT the thermal energy.
! When a defect has a total fixed concentration, the individual
! C_X^q (concentration of defect X in charge q) are calculated as:
! C_X^q = B_q * (C_X^all - sum_q' C_X^q') / (B_0 + sum_q/=q' B_q)
!
! Here, q' symbolises any charge states that are kept fixed.
! If C_X^0 is not been fixed, it is calculated as:
! C_X^0 = C_X^all - sum_q C_X^q - sum_q' C_X^q'
!
! If the total concentration is not fixed, then the non-fixed 
! concentrations are calculated as normal
! C_X^q = N_X B_q
! and the the total C_X^all determined by summing the individual
! contributions.
!
! J. Buckeridge June 2016, June 2017 for largefac
!
! UPDATE - May 2017
! BUGFIX - Feb 2018
! BUGFIX - Apr 2018
! BUGFIX - Apr 2018
! BUGFIX - Aug 2019
! BUGFIX - Jan 2020
! 
!
  implicit none
!
  integer, parameter :: nmax = 25000
  integer, parameter :: nless = 500
  integer, parameter :: bigint = 1e6
  real*8, parameter :: small = 1.d-14
  real*8, parameter :: lsmall = 1.d-7
  real*8, parameter :: vlarge = 1.014232054735d304
  real*8, parameter :: largefac = 700.d0
  real*8, parameter :: kboltz = 8.617343d-5 ! Boltzmann constant (eV/K)
  integer i, j, k, l, m, stat, stat1, stat2
  integer count, lownum, kcount, transnum(nless), linelength, pspace
  integer ndefects, nspinpol, ncharge(nless), nsite(nless), numdos
  integer kmax, kmin
  integer nfrozen, nfnew, countfroz(nless)
  integer iconc(nless), iconcch(nless)
  integer countdos
  integer highnum
  real*8 nelect, egap, charge(nless,nless), energy(nless,nless), deg(nless,nless)
  real*8 temperature, ener0(nless), deg0(nless), enerq
  real*8 pi, sum1, sum2, temp1, emin, emax, loweng, fermitrans(nless), enertmp
  real*8 lattvec(3,3), factor, volume
  real*8 edos(nmax), dos(nmax), dosu(nmax), dosd(nmax), de, expfac, expfac2
  real*8 delphi, estep, direction, n0, p0, lhs, rhs, diff, diffold
  real*8 rhstest, lhstest
  real*8 highener, sum3
  real*8 delphi_min
  real*8 conc0(nless), concd(nless,nless), concall(nless)
  real*8 concf(nless), concfn(nless), chargefn(nless)
  character*3 str1
  character*10 name(nless), namef(nless), namefn(nless)
  character*80 line
  logical check0(nless), check_under, check_under2, checkfroz, checkzero
  logical fixconc(nless), fixconcch(nless), fixcharg(nless,nless)
  logical checkconv
  logical checkposcar
  logical checkfixchg, checkfixchgdir
!
  pi = acos(-1.0)
!
! Write header
!
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') "  (FROZEN)"
  write(*,'(a)') 
  write(*,'(a)') "   SSSS    CCCC      FFFFFF  EEEEEE   RRRR   MM     MM  IIIII" 
  write(*,'(a)') "  SS   S  CC   C     FF      EE      RR   R  MMM   MMM    I" 
  write(*,'(a)') "  SS      CC         FF      EE      RR  R   M MM MM M    I" 
  write(*,'(a)') "   SSSS   CC     --- FFFFFF  EEEEEE  RRRR    M  MMM  M    I" 
  write(*,'(a)') "      SS  CC         FF      EE      R   R   M   M   M    I" 
  write(*,'(a)') "  S   SS  CC   C     FF      EE      R   RR  M   M   M    I" 
  write(*,'(a)') "   SSSS    CCCC      FF      EEEEEE  R   RR  M   M   M  IIIII" 
  write(*,'(a)') 
  write(*,'(a)') 
  write(*,'(a)') "Energies in eV, temperature in Kelvin, DOS in states/unit cell" 
  write(*,'(a)') 
  write(*,'(a)') "------"
  write(*,'(a)') "j.buckeridge@ucl.ac.uk 2018"
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
!
! Open unitcell.dat to read in vectors and determine volume of pure cell
!
  checkposcar = .false.
  open(unit=11,file="unitcell.dat",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "ERROR: unitcell.dat (unit cell for which DOS determined) is needed!!"
     checkposcar = .true.
  else
     write(*,'(a)') "unitcell.dat found..."
     write(*,'(a)') "(Should be cell for which DOS was determined!)"
     write(*,*)
!
! Read in scaling factor from the file
!
     count = 0
     do
        count = count + 1
        read(11,*) line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) factor
        if(stat /= 0) then
           write(*,'(a)') "ERROR: scaling factor for lattice vectors not found!!"
           stop
        else
           exit
        endif
        if(count > bigint) then
           write(*,'(a)') "ERROR: cannot read data in unitcell.dat!!"
           stop
        endif
     enddo
!
! Now read in the lattice vectors
!
     count = 0
     i = 0
     do
        count = count + 1
        read(11,*) line
        if(line(1:1) == "#") then
           cycle
        else
           backspace(unit=11)
           i = i+1
           read(11,*,iostat=stat) (lattvec(i,j),j=1,3)
           if(stat /= 0) then
              write(*,'(a39,x,i1,x,a2)') "ERROR: could not read in lattice vector", i, "!!"
              stop
           endif
           if(i == 3) exit
        endif
        if(count > bigint) then
           write(*,'(a)') "ERROR: lattice vectors not found in unitcell.dat!!"
           stop
        endif
     enddo
     lattvec = lattvec * factor
     close(unit=11)
  endif

!
! Open POSCAR to read in vectors and determine volume of pure cell
!
!
  if(checkposcar) then
     open(unit=11,file="POSCAR",status="old",iostat=stat)
!
     if(stat /= 0) then
        write(*,'(a)') "Exiting as no unit cell data can be found..."
        stop
     endif
     write(*,'(a)') "...but found a POSCAR! Will use that instead..."
     write(*,'(a)') "(Should be cell for which DOS was determined!)"
     write(*,*)
!
! Read in lattice vectors
!
     read(11,*)
     read(11,*) factor
     do i=1,3
        read(11,*) (lattvec(i,j),j=1,3)
     enddo
     lattvec = lattvec * factor
     close(unit=11)
  endif
!
! Calculate cell volume
!
  volume = lattvec(1,1) * (lattvec(2,2)*lattvec(3,3) - lattvec(2,3)*lattvec(3,2))&
       & + lattvec(1,2) * (lattvec(2,3)*lattvec(3,1) - lattvec(2,1)*lattvec(3,3))&
       & + lattvec(1,3) * (lattvec(2,1)*lattvec(3,2) - lattvec(2,2)*lattvec(3,1))
  write(*,'(a15,2x,f13.6,x,a3)') "Volume of cell:", volume, "A^3"
  write(*,*)
!
! Open input file and read in data
!
  open(unit=11, file="input-fermi-frozen.dat", status="old", iostat=stat)
  if(stat /= 0) then
     write(*,'(a)') "No input-fermi-frozen.dat file found - bye!!"
     stop
  endif
!
! Read in whether system is spin polarised or not (will be used 
! for reading in DOS)
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nspinpol
     if(stat /= 0) then
        write(*,'(a)') "ERROR: spin polarisation specification not found!!"
        stop
     else
        if(nspinpol == 1) then
           write(*,'(a)') "Found non-spin polarised system..."
           exit
        elseif(nspinpol == 2) then
           write(*,'(a)') "Found non-spin polarised system..."
           exit
        else
           write(*,'(a)') "ERROR: spin polarisation specification must by 1 or 2!!"
           stop
        endif
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: spin polarisation specification not found!!"
        stop
     endif
  enddo
  write(*,*)
!
! Read in number of electrons (will be used to renormalise DOS)
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nelect
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of electrons not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6)') "Number of electrons in system:", nelect
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of electrons not found!!"
        stop
     endif
  enddo
!
! Read in energy gap
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) egap
     if(stat /= 0) then
        write(*,'(a)') "ERROR: energy gap not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6,x,a2)') "Energy gap of system:", egap, "eV"
        if(egap < 0.d0) then
           write(*,'(a)') "You have a negative gap - this is going to get weird..."
        endif
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: energy gap not found!!"
        stop
     endif
  enddo
!
! Read in temperature
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) temperature
     if(stat /= 0) then
        write(*,'(a)') "ERROR: temperature not found!!"
        stop
     else
        write(*,'(a30,2x,f13.6,x,a1)') "Temperature :", temperature, "K"
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: temperature not found!!"
        stop
     endif
  enddo
!
! Read in number of defect species
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) ndefects
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of defect species not found!!"
        stop
     else
        write(*,'(a30,2x,i4)') "Number of defect species:", ndefects
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of defect species not found!!"
        stop
     endif
     if(temperature <= small) then
        write(*,'(a)') "ERROR: zero temperature means zero defects!!"
        stop
     endif
  enddo
  write(*,*)
!
! For each species read in name, number of sites in the unit cell, 
! charge states, degeneracy and energies
!
  if(ndefects > 0) then
     count = 0
     i = 1
     do
        count = count + 1
        read(11,'(a)') line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) name(i), ncharge(i), nsite(i)
        if(stat /= 0) then
           write(*,'(a)') "ERROR: failed to read in defect name/charge!!"
           stop
        else
           if(ncharge(i) <= 0) then
              write(*,'(a13,x,i4,x,a36)') "ERROR: defect", i, &
                   &"has idiotic number of charge states!!"
              stop
           endif
           j = 1
           do
              read(11,'(a)') line
              if(line(1:1) == "#") cycle
              read(line,*,iostat=stat) charge(i,j), energy(i,j), deg(i,j)
              if(stat /= 0) then
                 write(*,'(a29,i3,x,a4,i3,a2)') "ERROR: failed read for defect", &
                      &i, "line", j, "!!"
                 stop
              endif
              j = j+1
              if(j > ncharge(i)) exit
           enddo
        endif
        if(i >= ndefects) exit
        i = i+1
        if(count > bigint) then
           write(*,'(a41,i3,a2)') "ERROR: failed to find all data for defect", i, "!!"
           stop
        endif
     enddo
!
! For sanity check return the defects and charge states found
!
     write(*,'(a)') "Defects found:"
     write(*,'(a4,2x,a15,2x,a20)') "Name", "# charge states", "# sites in unit cell"
     do i=1,ndefects
        write(*,'(a10,7x,i4,18x,i4)') name(i), ncharge(i), nsite(i)
     enddo
     write(*,*)
  else
     write(*,'(a)') "0 defects found..."
     write(*,*)
  endif
!
! Read in frozen concentrations (in cm^-3). First read in number of
! defects with fixed overall charge
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nfrozen
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of fixed concentration defect species not found!!"
        stop
     else
        write(*,'(a45,2x,i4)') "Number of fixed concentration defect species:", nfrozen
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of fixed concentration defect species not found!!"
        stop
     endif
  enddo
!
! Initialise array that tracks which defects have fixed total concentrations
!
  fixconc = .false.
!
! For each of these species read in name and fixed total concentration, and 
! record which of the original defects they relate to
!
  if(nfrozen > 0) then
     count = 0
     countfroz = 0
     i = 1
     do
        count = count + 1
        read(11,'(a)') line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) namef(i), concf(i)
        if(stat /= 0) then
           write(*,'(a)') "ERROR: failed to read in frozen defect name/charge!!"
           stop
        else
!
! checkfroz is to see if frozen defect is part of the original list, iconc records
! which frozen defect relates to which original defect, countfroz is for checking
! whether a single defect has been assigned more than one frozen concentration
!
           checkfroz = .true.
           do j=1,ndefects
              if(name(j) == namef(i)) then
                 checkfroz = .false.
                 fixconc(j) = .true.
                 iconc(j) = i
                 countfroz(j) = countfroz(j) + 1
              endif
           enddo
           if(checkfroz) then
              write(*,'(a33,x,i3,x,a35)') "ERROR: fixed concentration defect", i, &
                   &"is not in list of original defects!"
              write(*,'(a)') "       please add new frozen defects as fixed charge defects in"
              write(*,'(a)') "       separate list - otherwise solution cannot be found!!"
              stop
           endif
        endif
        if(i >= nfrozen) exit
        i = i+1
        if(count > bigint) then
           write(*,'(a48,i3,a2)') "ERROR: failed to find all data for frozen defect", i, "!!"
           stop
        endif
     enddo
!
! Check if any defects have been assigned multiple fixed concentrations
! in error
!
     do i=1,ndefects
        if(countfroz(i) > 1) then
           write(*,'(a13,x,i3,x,a35)') "ERROR: defect", i, &
                &"has multiple fixed concentrations!!"
           write(*,'(a)') "       Nonsensical input - bye!"
           stop
        endif
     enddo
  endif
!
! Now read in number of species with fixed charged defect concentrations
!
  count = 0
  do
     count = count + 1
     read(11,*) line
     if(line(1:1) == "#") cycle
     read(line,*,iostat=stat) nfnew
     if(stat /= 0) then
        write(*,'(a)') "ERROR: number of fixed charged defect species not found!!"
        stop
     else
        write(*,'(a39,8x,i4)') "Number of fixed charged defect species:", nfnew
        exit
     endif
     if(count > bigint) then
        write(*,'(a)') "ERROR: number of fixed concentration defect species not found!!"
        stop
     endif
  enddo
  write(*,*)
!
! Initialise arrays fixconcch, which is for checking if a defect has a frozen
! charge state; and fixcharg, which tracks which defect charge state is frozen
!
  fixcharg = .false.
  fixconcch = .false.
!
! For each species read in name, charge and concentration. Check if they
! correspond to any defects in the original list
!
  checkfixchg = .false.
  if(nfnew > 0) then
     count = 0
     iconcch = 0
     i = 1
     do
        count = count + 1
        read(11,'(a)') line
        if(line(1:1) == "#") cycle
        read(line,*,iostat=stat) namefn(i), chargefn(i), concfn(i)
!
! If there is any fixed C_X^q with q/=0, we need to know for the small 
! exponentials routine later (see 19/4/18)
!
        if(abs(chargefn(i)) > small) checkfixchg = .true.
        if(stat /= 0) then
           write(*,'(a)') "ERROR: failed to read in fixed charge defect name/charge!!"
           stop
        else
           do j=1,ndefects
              if(name(j) == namefn(i)) then
!
! iconcch keeps track of which original defect each fixed charged defect relates to
!
                 fixconcch(j) = .true.
                 iconcch(i) = j
                 do k=1,ncharge(j)
                    if(abs(chargefn(i) - charge(j,k)) < small) fixcharg(j,k) = .true.
                 enddo
              endif
           enddo
        endif
        if(i >= nfnew) exit
        i = i+1
        if(count > bigint) then
           write(*,'(a48,i3,a2)') "ERROR: failed to find all data for frozen defect", i, "!!"
           stop
        endif
     enddo
  endif
!
! Now check for any defect which has fixed total concentration and a fixed
! concentration of a particular charge - if the concentration of the charge state(s)
! is greater than the total the input makes no sense. Also, if all charges have
! been fixed but do not sum up to the total, that also makes no sense
!
  do i=1,ndefects
     if(fixconcch(i)) then
        if(fixconc(i)) then
!
! First check for summed fixed charge states that is greater than the 
! fixed total charge
!
           sum1 = 0.d0
           do j=1,nfnew
              if(iconcch(j) == i) sum1 = sum1 + concfn(j)
           enddo
           if(sum1 - concf(i) > small) then
              write(*,'(a13,i3,x,a41)') "ERROR: defect", i, "has a fixed &
                   &total concentration less than"
              write(*,'(a)') "       that of its fixed concentration &
                   &of charge states!!"
              write(*,'(a)') "       Input is nonsensical, bye!!"
              stop
           endif
!
! Now check for cases where all charge states are fixed, but the sum is
! not equal to the total fixed concentration
!
           checkfroz = .true.
           do j=1,ncharge(i)
              if(fixcharg(i,j) .eqv. .false.) checkfroz = .false.
           enddo
           if(checkfroz) then
              sum1 = 0.d0
              do j=1,nfnew
                 if(iconcch(j) == i) sum1 = sum1 + concfn(j)
              enddo
              if(abs(sum1 - concf(i)) > small) then
                 write(*,'(a13,i3,x,a39)') "ERROR: defect", i, &
                      &"has all charge states fixed but the sum"
                 write(*,'(a)') "       of the concentrations &
                      &does not equal the total fixed"
                 write(*,'(a)') "       concentration!!"
                 write(*,'(a)') "       Input is nonsensical, bye!!"
                 stop
              endif
           endif
        endif
     endif
  enddo
!
! Write to screen the fixed defect concentrations as a check for the user
!
  if(nfrozen > 0 .or. nfnew > 0) then
     write(*,'(a)') "Frozen defects: (name)  (charge)   (concentration / cm^-3)"
     if(nfrozen > 0) then
        do i=1,nfrozen
           write(*,'(a16,x,a6,6x,a3,2x,e22.13e3)') "                ", &
                namef(i), "all", concf(i) 
        enddo
     endif
     if(nfnew > 0) then
        do i=1,nfnew
           write(*,'(a16,x,a6,f9.1,2x,e22.13e3)') "                ", &
                &namefn(i), chargefn(i), concfn(i)
        enddo
     endif
  endif
  write(*,*)
!
  close(unit=11)
!
! Read in density of states and renormalise. The energy scale should 
! be relative to the VBM (or Fermi level for a metal). Take into account
! spin polarisation specification
!
  open(unit=11, file="totdos.dat", status="old", iostat=stat)
  if(stat /= 0) then
     write(*,'(a)') "ERROR: no totdos.dat file found - bye!!"
     stop
  endif
!
  i = 0
  countdos = 0
  do
     read(11,'(a)',iostat=stat) line
     if(stat < 0) exit
     if(line(1:1) == "#") cycle
     i = i+1
     if(nspinpol == 1) then
        read(line,*) edos(i), dos(i)
        if(dos(i) < 0.d0) then
           countdos = countdos + 1
        endif
     else
        read(line,*) edos(i), dosu(i), dosd(i)
        if(dosu(i) < 0.d0 .or. dosd(i) < 0.d0) then
           countdos = countdos + 1
        endif
        dos(i) = dosu(i) + dosd(i)
     endif
  enddo
  close(unit=11)
  if(countdos > 0) then
     write(*,'(a)') "WARNING: Found negative value(s) of DOS!"
     write(*,'(a)') "         Do you know what you are doing?"
     write(*,'(a)') "         These may cause serious problems..."
     write(*,*)
  endif
  numdos = i
  emin = edos(1)
  emax = edos(numdos)
  if(egap > emax) then
     write(*,'(a)') "ERROR: your conduction band is not present in energy range of DOS!!"
     stop
  endif
!
! Now integrate DOS up to VBM (using trapezoidal rule) to get number of 
! states per unit cell. Renormalise according to number of electrons
!
  sum1 = 0.d0
  i=0
  do
     i = i+1
     if(edos(i+1) > 0.d0) exit
     de = (edos(i+1) - edos(i))
     sum1 = sum1  + (dos(i) + dos(i+1)) * de / 2.d0
     if(i+1 > numdos) then
        write(*,'(a)') "ERROR: Energy scale in DOS does not cross zero!!"
        stop
     endif
  enddo
  write(*,'(a37,2x,f13.6)') "Integration of DOS up to Fermi level:", sum1
  dos = dos * nelect / sum1
  sum1 = 0.d0
  i=0
  do
     i = i+1
     if(edos(i+1) > 0.d0) exit
     de = (edos(i+1) - edos(i))
     sum1 = sum1 + (dos(i) + dos(i+1)) * de / 2.d0
     if(i+1 > numdos) then
        write(*,'(a)') "ERROR: Energy scale in DOS does not cross zero!!"
        stop
     endif
  enddo
  write(*,'(a37,2x,f13.6)') "Renormalised integrated DOS         :", sum1
  write(*,*)
!
! Now do algorithm to determine self-consistent Fermi level. Begin scanning
! at E_F = 0, go forwards by 1/10 of difference between E_F = 0 and E_F = 
! emax. Determine difference of (*) from zero, when sign of difference changes,
! reverse direction of scan reducing step size. Continue until step
! size is lower than a particular convergence criterion, or the difference
! is less than same. Take care for large exponentials (which return infinity),
! if they occur set the calculated value at a max. If this max is maintained
! between more than two steps, reverse the direction of the scan.
!
!
! Calculate neutral defect concentrations. Skip if the defect has fixed total
! concentration, as the neutral concentration must be calculated aferwards
! in that case. Also check for fixed neutral states.
!
  if(ndefects > 0) then
     conc0 = 0.d0
     do i=1,ndefects
        check0(i) = .true.
        do j=1,ncharge(i)
           if(abs(charge(i,j)) < small) then
              ener0(i) = energy(i,j)
              deg0(i) = deg(i,j)
              if(fixcharg(i,j)) then
                 check0(i) = .true.
              else
                 check0(i) = .false.
              endif
           endif
        enddo
        if(check0(i) .or. fixconc(i)) cycle
        expfac = -ener0(i) / (kboltz * temperature)
        if(expfac > largefac) then
           conc0(i) = vlarge
        else
           conc0(i) = nsite(i) * deg0(i) * exp(expfac)
        endif
     enddo
  endif
!
! Set initial values for energy scan
!
  estep = emax / 10.d0
  direction = 1.d0
  delphi = -estep
  check_under = .false.
  check_under2 = .false.
  concd = 0.d0
!
  i = 0
  kmax = 0
  kmin = 0
  checkconv = .false.
  checkfixchgdir = .false.
!
  do
!
! Step E_F (don't step E_F if convergence in estep has been found in the
! previous iteration)
!
     if(checkconv .eqv. .false.) delphi = delphi + direction * estep
!
! Get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
!
     n0 = 0.d0
     p0 = 0.d0
     do j=1,numdos-1
        if(edos(j+1) <= 0.d0) then
           de = (edos(j+1) - edos(j))
           p0 = p0 + ( dos(j) / ( 1 + exp((delphi - edos(j)) / (kboltz * temperature)) )&
                & + dos(j+1) / (1 + exp((delphi - edos(j+1)) / (kboltz * temperature)) ) )&
                * de / 2.d0
        endif
        if(edos(j) >= egap) then
           de = (edos(j+1) - edos(j))
           n0 = n0 + ( dos(j) / ( 1 + exp((edos(j) - delphi) / (kboltz * temperature)) )&
                & + dos(j+1) / (1 + exp((edos(j+1) - delphi) / (kboltz * temperature)) ) )&
                * de / 2.d0
        endif
     enddo
     lhs = p0
     rhs = n0
!
! Get defect concentrations at E_F. If the total concentration is fixed,
! then get C_X^q first (taking into account any charge states that are fixed),
! if it isn't fixed then use normal formula
!
     if(ndefects > 0) then
        do j=1,ndefects
           if(fixconc(j)) then
!
! Sum up Gibbs factors B_q = g_q exp(-E_f^q/kT), skipping q=0 and any fixed q'
! BUGFIX 17/04/18****
! Note for low T and high E_f^q, exponentials can end up as zero, and the
! distribution of the total fixed concentration amongst the various charge states
! becomes impossible. Instead, we must find the lowest E_f^q at the current delphi
!
              loweng = energy(j,1) + charge(j,1) * delphi
              lownum = 1
              if(ncharge(j) > 1) then
                 do k=1,ncharge(j)
                    enertmp = energy(j,k) + charge(j,k) * delphi
                    if(enertmp - loweng <= small) then
                       loweng = enertmp
                       lownum = k
                    endif
                 enddo
              endif
!
! Now we take B_q = g_q / g_qmin exp(-(E_f^q - E_f^qmin)/kT, meaning we divide
! above and below by the B_q factor for the lowest formation energy state. We
! thus avoid zero summations as there will always be a term equal to one (unless
! the T is ridiculously low, and in that case we set everything to zero to avoid
! infinities). We get the sum_q/=0 B_q here:
!
              sum1 = 0.d0
              do k=1,ncharge(j)
                 if(fixcharg(j,k) .or. abs(charge(j,k)) < small) cycle
                 enerq = energy(j,k) + charge(j,k) * delphi - loweng
                 expfac = -enerq / (kboltz * temperature)
                 if(expfac > largefac) then
                    sum1 = sum1 + vlarge
                 else
                    sum1 = sum1 + (deg(j,k) / deg(j,lownum)) * exp(expfac)
                 endif
              enddo
!
! Get fixed C_X^all and subtract all C_X^q' (where q' is a fixed charge,
! including q=0)
!
              sum2 = concf(iconc(j)) * volume / 1.d24
              if(fixconcch(j)) then
                 do k=1,nfnew
                    if(iconcch(k) == j) then
                       sum2 = sum2 - concfn(k) * volume / 1.d24
                    endif
                 enddo
              endif
!
! BUGFIX 18/4/18****
! Now we check if the exponentials are still zero, despite using E_f^qmin
! which will be the case and problematic if sum1 is zero. If they are, we
! find the max energy, and distribute the fixed concentration among the 
! charge states (that can vary) according to:
! C_X^q = [E_f^qmax - E_f^q]^8 / sum_q/=q' [E_f^qmax - E_f^q]^8
!
              if(sum1 == 0.d0) then
!
! Find max energy at current delphi used for distribution
!
                 highener = energy(j,1) + charge(j,1) * delphi
                 highnum = 1
                 if(ncharge(j) > 1) then
                    do l=1,ncharge(j)
                       if(fixcharg(j,l)) cycle
                       enertmp = energy(j,l) + charge(j,l) * delphi
                       if(enertmp - highener > small) then
                          highener = enertmp
                          highnum = l
                       endif
                    enddo
                 endif
!
! BUGFIX 19/4/18****
! Add a small 'delta' to the max energy, so that the q with highest energy
! doesn't automatically get set a zero concentration (the zero concentration
! could cause problems when there are few charge states)
!
                 highener = highener + lsmall
              endif
!
! Now calculate C_X^q = B_q (C_X^all - sum_q' C_X^q') / (B_0 + sum_q/=q' B_q). 
! If the neutral state was not supplied or if neutral is fixed, then set 
! B_0 = 0 (as C_X^0 has been subtracted from C_X^all, this is correct for 
! the fixed neutral case. Where neutral was not supplied, B_0 = 0 means 
! E_f^0 is infinity)
!
              do k=1,ncharge(j)
                 if(fixcharg(j,k) .or. abs(charge(j,k)) < small) cycle
                 if(sum1 == 0.d0) then
!
! If exponentials are too small distribute C_X^q as described two comments above
!
                    sum3 = 0.d0
                    do l=1,ncharge(j)
                       if(fixcharg(j,l)) cycle
                       sum3 = sum3 + (highener - energy(j,l) - charge(j,l) * delphi)**8
                    enddo
                    concd(j,k) = (highener - energy(j,k) - charge(j,k) * delphi)**8 &
                         & / sum3
                    concd(j,k) = concd(j,k) * sum2
                 elseif(sum1 == 1.d0) then
!
! BUGFIX 24/4/18****
! If X^qmin dominates, so that the sum of B_q (sum1 above) is exactly equal
! to 1 (due to exp(-(E_f(X^q) - E_f(X^qmin))/kT) = 1), we need to set
! C_X^qmin = C_X^T - sum_q' C_X^q', and all other C_X^q = 0. Doing so avoids
! situation where other C_X^q are nonzero, but have not contributed to sum1
! because they are too small
!
! BUGFIX 14/1/20****
! Had problem with the fix above (24/4/18). When one charge dominated such that
! sum1 = 1.0, the 'fix' ignored the concentration of the neutral state (in cases
! where the neutral state was not specifically fixed, so needed to be evaluated
! as part of the total fixed concentration). So to correct for this, the B_0
! term has been added to the denominator when calculating the concentrations 
! of the individual charge states, i.e. we set
! C_X^qmin = (C_X^T - sum_q' C_X^q') / (B_0 + 1.0); where no q'=0
! All other C_X^q are still set to zero, as necessary (here q /= 0)
!
                    if(k == lownum) then
                       if(check0(j)) then
                          concd(j,k) = sum2
                       else
                          expfac2 = -(ener0(j) - loweng) / (kboltz * temperature)
                          if(expfac2 > largefac) then
                             concd(j,k) = sum2 / (vlarge + sum1)
                          else
                             concd(j,k) = sum2 / &
                                  &(deg0(j) * exp(expfac2) / deg(j,lownum) + sum1)
                          endif
                       endif
                    else
                       concd(j,k) = 0.d0
                    endif
                 else
                    enerq = energy(j,k) + charge(j,k) * delphi - loweng
                    expfac = -enerq / (kboltz * temperature)
                    if(check0(j)) then
                       if(expfac > largefac) then
                          concd(j,k) = vlarge * sum2 
                       else
                          concd(j,k) = (deg(j,k) / deg(j,lownum)) * exp(expfac) * sum2 
                       endif
!
! Now divide by sum1, i.e. sum_q/=q' B_q
!
                       concd(j,k) = concd(j,k) / sum1
!
                    else
                       expfac2 = -(ener0(j) - loweng) / (kboltz * temperature)
                       if(expfac > largefac) then
                          concd(j,k) = vlarge * sum2 
                       else
                          concd(j,k) = (deg(j,k) / deg(j,lownum)) * exp(expfac) * sum2 
                       endif
!
! Now divide by B_0 + sum_q/=q' B_q
!
                       concd(j,k) = concd(j,k) / &
                            &(deg0(j) * exp(expfac2) / deg(j,lownum) + sum1)
!
                    endif
                 endif
              enddo
           else
!
! Total concentration not fixed - get C_X^q = N_X B_q. C_X^0 calculated previously
! and zero if E_f^0 not supplied. Fixed charge states do not affect this part
!
              do k=1,ncharge(j)
                 if(fixcharg(j,k) .or. abs(charge(j,k)) < small) cycle 
                 enerq = energy(j,k) + charge(j,k) * delphi
                 expfac = -enerq / (kboltz * temperature)
                 if(expfac > largefac) then
                    concd(j,k) = vlarge
                 else
                    concd(j,k) = nsite(j) * deg(j,k) * exp(expfac)
                 endif
              enddo
           endif
!
! Now add all C_X^q to RHS or LHS as appropriate
!
           do k=1,ncharge(j)
              if(fixcharg(j,k)) cycle
              if(charge(j,k) < small) then
                 rhs = rhs + concd(j,k) * abs(charge(j,k))
              else
                 lhs = lhs + concd(j,k) * abs(charge(j,k))
              endif
           enddo
        enddo
     endif
!
! Before summing all concentrations to LHS and RHS, first we need to sum all 
! concentrations excluding those that have been fixed, in order to check 
! whether there is underflow - so that we can enter the routine that solves
! the case where the exponentials are set to zero and nothing changes in the 
! normal routine. In SC-FERMI, this check is done by seeing if LHS = RHS = 0
! (absolutely), but here we need to account for the possibility of the LHS
! and RHS not being zero, due to fixed concentrations, but still have nothing 
! changing as E_F changes (see 10/4/18). We exclude from this procedure those
! defects that have C_X^T fixed, as they can still have subdivided C_X^q that
! change (see 18/4/18)
!
! We set new variables LHStest and RHStest, so that the new test is if
! LHStest = RHStest = 0 (absolutely)
!
     rhstest = rhs
     lhstest = lhs
!
! Add on all fixed charge defects, whether part of the original set or not
! (it makes no difference - the important point is to change LHS or RHS).
! Sum them up first as you add to the RHS, i.e. positively summed if negatively
! charged and vice verca. The sum can then be added to the RHS afterwards.
! The point of this is to avoid having v. large numbers on RHS and LHS that
! exactly cancel, but make it practically impossible to find solution by 
! balancing the other relatively tiny numbers
!
     sum1 = 0.d0
     if(nfnew > 0) then
        do j=1,nfnew
           if(chargefn(j) < small) then
              sum1 = sum1 + abs(chargefn(j)) * concfn(j) * volume / 1.d24
           else
              sum1 = sum1 - abs(chargefn(j)) * concfn(j) * volume / 1.d24
           endif
        enddo
     endif
     rhs = rhs + sum1
!
! Get diff, check how small it is, get direction to go to reduce it
!
     diff = rhs - lhs
!
! If estep was converged in the previous iteration, exit now after recalculating
! all concentrations and diff (see 19/4/18)
!
     if(checkconv) exit
!
! Need to check for real underflow in exponentials, which is a problem if
! both the rhs and lhs are identically zero. If all terms become
! zero due to underflow, we find the range of fermi level within which
! this happens, and as an approximation return the midpoint as the solution
!
     if(rhstest == 0.d0 .and. lhstest == 0.d0) then
        if(check_under) cycle
!
! For the case where some C_X^q fixed and q/=0, we may need to keep scanning
! beyond region where everything stops changing due to small exponents. After
! coming back from such a continued scan, we need to go back to delphi_min to
! find boundaries of range (because solution must be in the range).
! Remember, because we have changed direction already, we continue in the same
! direction, and step delphi back one as that is the step just before everything
! went to zero (see 19/4/18)
!
        if(checkfixchgdir) then
           check_under = .true.
           check_under2 = .true.
           temp1 = delphi - direction * estep
           delphi = delphi_min
           delphi_min = temp1
           cycle
        endif
        check_under = .true.
        write(*,'(a)') "WARNING: real number underflow due to v. small exponentials"
        write(*,'(a)') "         - treat results with caution!!"
        delphi_min = delphi
        if(estep * 1.d-4 > small) then
           estep = estep * 1.d-4
        else
           estep = small
        endif
     else
        if(check_under) then
           if(check_under2) then
!
! BUGFIX 23/4/18****
! After getting midpoint of range where everything went to zero as solution,
! need to go back and recalculate concentrations and diff. Otherwise, in very
! rare cases can get concentrations that are too large (= vlarge, but when 
! multiplied by 10^24 to get correct units, end up with 'infinity'), as the
! concentrations were calculated after exponentials were no longer zero
!
              delphi = 0.5d0 * (delphi - direction * estep + delphi_min)
              checkconv = .true.
              cycle
           else
!
! BUGFIX 19/4/18****
! If some C_X^q q/=0 are fixed, we need to check if the solution is still in
! the forward direction, so instead of flipping back to delphi_min and getting
! the boundaries of the region where everything went to zero, we first keep 
! scanning (with the larger step size used before entering the region)
! see 19/4/18
!
              if(checkfixchg) then
                 if(check_under) estep = estep * 1.d4
                 check_under = .false.
                 checkfixchgdir = .true.
              else
                 check_under2 = .true.
                 temp1 = delphi
                 delphi = delphi_min - direction * estep
                 delphi_min = temp1
                 direction = direction * (-1.d0)
              endif
           endif
        else
!
! Otherwise all is ok and we can perform the normal procedure
!
! BUGFIX APR 2018: put an additional check for check_under2, to 
! avoid the situation where you get a direction change directly 
! after a direction change when check_under2 is set to true
! (18/4/18 + 19/4/18)
!
           if(abs(diff) <= small*(rhs + lhs)) exit
!
! BUGFIX SEP 2018 (see 14-9-18): removed block of code that was originally
! put in place to check if diff = diffold and diff >= vlarge. The method was
! cumbersome, relying on many checks, and failed sometimes. Now the code 
! deals with vlarge diffs in the following way:
! 1. If |diff| > vlarge but |diffold| < vlarge, delphi goes back to previous
!    point and steps in same direction as before with smaller estep (to 
!    find where diff gets too big along a finer grid and change direction
!    appropriately)
! 2. If |diff| > |diffold| but both are > vlarge, direction is changed and
!    step size is reduced, even if it occurs in the first two steps (to avoid
!    going off into a region where diff > vlarge and getting stuck)
! 3. If |diff| < |diffold|, both are > vlarge, but there is a sign change
!    between them (checked by dividing one by the other and seeing if its < 0)
!    the direction is changed and step size reduced as the solution lies 
!    somewhere between them (the same situation but with |diff| > |diffold| is
!    covered in 2 above)
!
           if(i > 0) then 
              if(abs(diff) > abs(diffold)) then
                 if(abs(diff) >= vlarge .and. abs(diffold) < vlarge) then
                    delphi = delphi - direction * estep
                    diff = diffold
                    estep = estep / 4.d0
                 else
                    direction = direction * (-1.d0)
                    if(i > 1 .or. abs(diff) >= vlarge) estep = estep / 10.d0
                 endif
              else
                 if( (abs(diff) > vlarge .and. abs(diffold) > vlarge) .and. &
                      (diff / diffold < 0.d0) ) then
                    direction = direction * (-1.d0)
                    estep = estep / 10.d0
                 endif
              endif
           endif
        endif
        if(estep < small) then
!
! The solution is deemed converged. However, we must take a step back and
! recalculate the concentrations at the previous delphi, as these will be
! closer to the true solution (diff is lower in the previous step). Remember
! that direction has been changed, so we add rather than subtract direction *
! estep * 10.d0 (the 10 accounts for the decrease in step size between steps)
! See 19/4/18
!
           checkconv = .true.
           delphi = delphi + direction * estep * 10.d0 
        endif
        i = i+1
        diffold = diff
!
! If scan gets to boundaries of energy scale reverse direction, but only once
!
        if(delphi > emax) then
           if(kmax == 0) then
              direction = direction * (-1.d0)
              kmax = kmax+1
           else
              write(*,'(a)') "Sorry, no solution found within range of DOS :("
              stop
           endif
        endif
        if(delphi < emin) then
           if(kmin == 0) then
              direction = direction * (-1.d0)
              kmin = kmin+1
           else
              write(*,'(a)') "Sorry, no solution found within range of DOS :("
              stop
           endif
        endif
     endif
  enddo
!
! Write results - E_F and concentrations
!
  write(*,'(a)') "Solution found!"
  write(*,*)
  write(*,'(a41,x,e22.13e3,2x,a5)') "Condition (n + acceptors) - (p + donors):", &
       &diff * 1.d24 / volume, "cm^-3"
  write(*,*) 
  write(*,'(a)') "Results:"
  write(*,'(a)') "--------"
  write(*,'(a16,x,f21.13,2x,a4)') "SC Fermi level :", delphi, "(eV)"
  write(*,*)  
  write(*,'(a)') "Concentrations:"
  write(*,'(a16,x,e22.13e3,2x,a5)') "n (electrons)  :", n0 * 1.d24 / volume, "cm^-3" 
  write(*,'(a16,x,e22.13e3,2x,a5)') "p (holes)      :", p0 * 1.d24 / volume, "cm^-3"
!
! Now for defects, make sure the fixed concentrations are accounted for
!
  if(ndefects > 0) then
     do i=1,ndefects
        if(fixconc(i)) then
!
! With fixed total concentration - set concall to that fixed concentration
!
           sum1 = concf(iconc(i)) * volume / 1.d24
           concall(i) = sum1
!
! Now need to get C_X^0 = C_X^all - sum_q C_X^q - sum_q' C_X^q' (where q' are
! fixed charge states). If q' = 0, then C_X^0 is a fixed concentration
!
           checkfroz = .true.
           if(fixconcch(i)) then
!
! Need to get sum_q' C_X^q', but check for q' = 0 (if so, set C_X^0 to the 
! fixed neutral concentration)
! 
              do j=1,nfnew
                 if(iconcch(j) == i) then
                    if(abs(chargefn(j)) < small) then
                       conc0(i) = concfn(j) * volume / 1.d24
                       checkfroz = .false.
                    else
                       sum1 = sum1 - concfn(j) * volume / 1.d24
                    endif
                 endif
              enddo
           endif
!
! Now need sum_q C_X^q (unless we had q' = 0). For any fixed q' concd = 0
!
           if(checkfroz) then
              sum2 = 0.d0
              do j=1,ncharge(i)
                 sum2 = sum2 + concd(i,j)
              enddo
!
! Get C_X^0. If neutral was not supplied originally, then sum1 = sum2. There may
! be slight errors from subtracting large numbers, so set to zero if the difference
! is very small
!
              conc0(i) = sum1 - sum2
              if(abs(conc0(i) / concall(i)) < small) conc0(i) = 0.d0
           endif
        else
!
! No fixed total concentration. Get the total by first getting C_X^0 (which will
! be zero if the neutral state was not supplied, or if it is a fixed charge state)
!
           concall(i) = conc0(i)
!
! Now add on the C_X^q
!
           do j=1,ncharge(i)
              concall(i) = concall(i) + concd(i,j)
           enddo
           if(fixconcch(i)) then
!
! If there are q', then add on the C_X^q'. If the neutral is fixed (i.e. q' = 0)
! then it is included here to C_X^all, but record its concentration in conc0
!
              sum1 = 0.d0
              do j=1,nfnew
                 if(iconcch(j) == i) then
                    sum1 = sum1 + concfn(j) * volume / 1.d24
                    if(abs(chargefn(j)) < small) conc0(i) = concfn(j) * volume / 1.d24
                 endif
              enddo
              concall(i) = concall(i) + sum1
           endif
        endif
!
! Write the total concentrations 
!
        write(*,'(a10,a6,x,e22.13e3,2x,a5)') name(i), "     :", concall(i) * 1.d24 / &
             &volume, "cm^-3"
     enddo
     write(*,*)
!
! Now the charge state breakdown. If C_X^all = 0 we can't do this part
!
     write(*,'(a)') "Breakdown of concentrations for each defect charge state:"
     do i=1,ndefects
        write(*,'(a)') "--------------------------------------------------------"
        if(concall(i) == 0.d0) then
           write(*,'(a10,a6,3x,a34)') name(i), "     :", "Zero total - &
                &cannot give breakdown"
           cycle
        endif
        write(*,'(a10,a6,3x,a6,3x,a20,x,a7)') name(i), "     :", "Charge", &
             &"Concentration(cm^-3)", "% total"
        do j=1,ncharge(i)
!
! Write the charge states that weren't fixed first
!
           if(fixcharg(i,j)) cycle
           if(abs(charge(i,j)) < small) then
              write(*,'(15x,a1,3x,f5.1,3x,e22.13e3,x,f7.2)') ":", charge(i,j), conc0(i) * &
                   1.d24 / volume, conc0(i) * 100.d0 / concall(i)
           else
              write(*,'(15x,a1,3x,f5.1,3x,e22.13e3,x,f7.2)') ":", charge(i,j), concd(i,j) * &
                   1.d24 / volume, concd(i,j) * 100.d0 / concall(i)
           endif
        enddo
        if(fixconcch(i)) then
!
! Now write the fixed charge states
!
           do j=1,nfnew
              if(iconcch(j) == i) then
                 write(*,'(15x,a1,3x,f5.1,3x,e22.13e3,x,f7.2)') ":", chargefn(j), &
                      &concfn(j), (concfn(j) * volume / 1.d24) * 100.d0 / concall(i)
              endif
           enddo
        endif
     enddo
  endif
  if(nfnew > 0) then
!
! Check for any fixed charged defects that weren't part of the original set
! and write their concentrations here
!
     count = 0
     do i=1,nfnew
        if(iconcch(i) == 0) then
           count = count + 1
           if(count == 1) then
              write(*,*)
              write(*,'(a)') "Concentrations of new fixed charge defects:"
              write(*,'(a10,a6,x,e22.13e3,2x,a5)') namefn(i), "     :", concfn(i), "cm^-3"
           else
              write(*,'(a10,a6,x,e22.13e3,2x,a5)') namefn(i), "     :", concfn(i), "cm^-3"
           endif
        endif
     enddo
  endif
!
end program frozen_fermi
