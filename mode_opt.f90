module GeoOptNormalCoordsApplicationModule

   !This module can be used for optimizations in normal coordinates.
   ![1] P. Bour and T.A. Keiderling, J. Chem. Phys. 117, 4126 (2002); doi: 10.1063/1.1498468

#include <DontWarnUnused.fh>

   use AMSApplicationModule
   use AMSResultsModule
   use AMSSolveRequestModule
   use ftlStringModule
   use InputReaderModule
   use KF
   use KFHistoryModule
   use NewNormalModesModule
   use VarTypes

   implicit none
   private

   type, extends(AMSApplicationType), public :: GeoOptNormalCoordsApplicationType
      integer(KINT)     :: maxIter
      integer(KINT)     :: iStep
      logical           :: calcHessian
      real(KREAL)       :: gradConv
      real(KREAL)       :: highFreq
      real(KREAL)       :: lowFreq
      type(ftlString)   :: hessSet
      type(ftlString)   :: fileName
   contains
      private
      ! Implementation of the AMSApplication interface
      procedure, public :: NewDerived
      procedure, public :: PrintHeaderDerived
      procedure, public :: RunDerived
      procedure, public :: DeleteDerived
   end type

   real(KREAL), parameter :: PI = 4.0_KREAL * atan(1.0_KREAL)

contains


   subroutine NewDerived(self, input)

      class(GeoOptNormalCoordsApplicationType)  , intent(inout) :: self
      type(InputType)                     , intent(inout) :: input

      
      call input%Get('GeoOptNormalCoords%MaxFreq'     , self%highFreq)    !Upper bound of frequency range
      call input%Get('GeoOptNormalCoords%MinFreq'     , self%lowFreq)     !Lower bound of frequency range
      call input%Get('GeoOptNormalCoords%MaxIter'     , self%maxIter)     !Maximum number of iterations
      call input%Get('GeoOptNormalCoords%NMGrad'      , self%gradConv)    !Convergence threshold for the gradients in normal coords
      call input%Get('GeoOptNormalCoords%InitHess'    , self%hessSet, convertCase="L") !Choose how to obtain the initial Hessian
      if (self%hessSet == 'file') then   !Retrieves Hessian from RKF
         if (.not. input%IsPresent('GeoOptNormalCoords%HessPath')) then
            call Stopit('ERROR: AMS: GeoOptNormalCoords: User did not supply initial Hessian.')
         endif
         call input%Get('GeoOptNormalCoords%HessPath'       , self%fileName)    !Leads to rkf file in path provided by user
      endif
      call input%Get('GeoOptNormalCoords%CalcHessUpdate'    , self%calcHessian) !Whether to calculate the Hessian at each iteration

   end subroutine



   subroutine PrintHeaderDerived(self, unit)

      class(GeoOptNormalCoordsApplicationType)  , intent(in) :: self
      integer(KINT)                       , intent(in) :: unit

      write(unit, '(A)') '       ***********************************************'
      write(unit, '(A)') '       * GEOMETRY OPTIMIZATION IN NORMAL COORDINATES *'
      write(unit, '(A)') '       ***********************************************'

      call self%chemicalSystem%Print(iuout)
      write(unit, '(A,t50,f10.2,A)') 'Lower bound of the frequency range:'             , self%lowFreq,  ' cm^-1'
      write(unit, '(A,t50,f10.2,A)') 'Higher bound of the frequency range:'            , self%highFreq, ' cm^-1'
      write(unit, '(A,t52,e8.2,A)') 'Convergence threshold for normal mode gradients:' , self%gradConv, ' Hartree/Bohr'
      write(unit, '(A,t51,I4)') 'Maximum number of iterations:'                        , self%maxIter
      write(unit, '(A,t52,L)') 'Calculate Hessian at every step?:'                     , self%calcHessian
      write(unit, '(A,t52,A)') 'Initial Hessian from:'                                 , self%hessSet%raw
      if (self%hessSet == 'file') write(unit, '(A)') 'Path for Hessian:'//self%fileName

   end subroutine



   subroutine RunDerived(self)

      class(GeoOptNormalCoordsApplicationType), intent(inout) :: self
      class(AMSResultsType)         , allocatable     :: results, prevResults
      integer(KINT)                                   :: iuKF
      logical                                         :: converged
      logical                       , external        :: master
      real(KREAL)                   , allocatable     :: cartHessian(:,:), NMGradients(:,:), modesInRange(:,:)
      real(KREAL)                   , external        :: convrs
      type(AMSSolveRequestType)                       :: request
      type(KFHistoryType)                             :: kfHistory
      type(NormalModesType)                           :: newModes


      call msg('Optimization in normal coordinates has started.')

      !Converts input frequency range to Hartree
      self%lowFreq  = self%lowFreq / (2.0_KREAL * PI * convrs('ENERGY','HARTREE,CM1'))
      self%highFreq = self%highFreq / (2.0_KREAL * PI * convrs('ENERGY','HARTREE,CM1'))

      !The master process needs to set up some things on the ams.rkf file, so that the GUI can visualize our results.
      if (master()) then
         call kfwrite(self%kfUnit, 'General%termination status', 'IN PROGRESS')  !For visualization while the job is running.
         call kfHistory%New(self%kfUnit, 'History')                              !Opens the History section from which the GUI reads.
         call kfwrfl(self%kfUnit)                                                !Flushes these changes to disk.
      end if

      request%title        = 'Initial'     !First calculation is carried out at given geometry to obtain the normal coordinates.
      request%gradients    = .true.        !Gradients are needed from now on.
      request%quiet        = .true.        !Asks the engine not to print to stdout, we get way too much output otherwise.
      if (self%hessSet%raw == 'calculate') request%normalModes = .true.  !Asks the engine to calculate normal modes.

      !Invokes the engine, passing the current geometry and the request. This will allocate and populate the results object.
      !The prevResults is optional in the engine (unallocated counts as not present). If it is present the engine can use the
      !old results to restart for example the SCF cycle, or whatever this specific engine needs to do internally.
      call self%engine%Solve(self%chemicalSystem, request, results, prevResults)
   
      !Checks if anything went wrong on the engine.
      if (.not.results%success) call StopIt('Engine failed to solve at initial geometry.')

      !Whether to calculate the Hessian at every step (True) or use BFGSupdate (False)
      request%normalModes = self%calcHessian 

      if (self%hessSet%raw == 'calculate') then
         !Calculated modes within the selected range are obtained as modesInRange
         call modeSetup(self, results%normalModes, modesInRange)
         !Saves Cartesian Hessian with the same name as above for simplification  
         cartHessian = results%hessian
      else
         if (self%hessSet == 'file') then
            !Load Hessian from RKF
            call loadHessian(self, cartHessian)
         else
            !Obtains a guess Hessian in Cartesian coordinates
            call initHessian(self, cartHessian)
         end if
         !After generating the guess Hessian or getting it from file, the initial modes are set up
         call newModes%New(self%chemicalSystem, cartHessian)
         !Modes within the selected range are then obtained as modesInRange
         call modeSetup(self, newModes, modesInRange)
      endif

      !Number of modes initially found in the range is printed to the output (this may change as the Hessian is updated)
      write (iuout,'(A,I5)') 'Number of normal modes (initially) found in the selected range', size(modesInRange(1,:))
      write (iuout,'(A)')

      !Calculate normal mode gradients for the first step in the optimization
      !Eqn. 6 (Step 4) in Reference [1]
      !NMGradients: g = S^t \cdot g_c; [nModes x 1] = [nModes x 3*nAtoms] * [3*nAtoms x 1]
      NMGradients = matmul(transpose(modesInRange), reshape(results%gradients,([3*self%chemicalSystem%nAtoms,1])))

      !Master writes the geometry along with the energy, gradients, etc. to the ams.rkf file, so that it can be viewed in the GUI.
      if (master()) then
         call kfHistory%Write('Coords', self%chemicalSystem%xyzAtoms)
         call kfHistory%Write('Energy', results%energy)
         call kfHistory%Write('Gradients', results%gradients)
         call kfHistory%Write('maxGrad', maxval(abs(results%gradients)))
         call kfHistory%Write('rmsGrad', sqrt( sum(results%gradients**2) / real(size(results%gradients),KREAL) ))
         call kfHistory%Write('nSelecNormalModes', size(modesInRange(1,:)))
         call kfHistory%Write('NMGradients', NMGradients)
         call kfHistory%FinishEntry()
      endif

      !Prints initial data in the output.
      write (iuout,'(A,f16.8,2(A,f12.8))') 'Step:    0, Energy = ', results%energy, ', MaxNMGrad = ', &
                                             maxval(abs(NMGradients)), ', MaxGrad = ', maxval(abs(results%gradients))

      !Check convergence and adjust the geometry if we are not converged yet.
      converged = all(abs(NMGradients) < self%gradConv)

      if (converged) then
         write (iuout, '(A)') 'Geometry already optimized.'
         deallocate(modesInRange)
         deallocate(NMGradients) 
      else         
         !Move results to prevResults
         call move_alloc(results, prevResults)

         !Displace the structure along the mode for the first step
         call newStep(self, cartHessian, NMGradients, modesInRange)

         !Deallocates as the size of it may change during optimization, will be allocated at every step
         deallocate(modesInRange)
         deallocate(NMGradients) 

         !Starts the loop with the new geometry generated in the previous step
         self%iStep = 0
         !Iterates until converged
         do while (.not. converged)

            self%iStep = self%iStep + 1

            !Breaks the loop if maxIter+1 is achieved
            if (self%iStep > self%maxIter) then
               write (iuout, '(A)') 'Number of iterations reached limited before convergence.'
               exit 
            endif

            !Gives every invocation of the engine a unique name.
            request%title = 'Step'//ftlString(self%iStep)

            !Calls the engine again with the new geometry and checks run.
            call self%engine%Solve(self%chemicalSystem, request, results, prevResults)
            if (.not.results%success) call StopIt('Engine failed to solve at step '//ftlString(self%iStep)//'.')
            
            !Obtains new Cartesian Hessian
            if (.not. self%calcHessian) then
               !Broyden–Fletcher–Goldfarb–Shanno Hessian update algorithm is applied to generate new Cartesian Hessian
               call BFGSupdate(self, prevResults, results, cartHessian)
               !Gets modes
               call newModes%New(self%chemicalSystem, cartHessian)
               !Gets modes within the selected range
               call modeSetup(self, newModes, modesInRange)
            else if (self%calcHessian) then 
               !Gets modes within the selected range
               call modeSetup(self, results%normalModes, modesInRange)
               !Save Cartesian Hessian in results%hessian for consistency with the guess Hessian     
               cartHessian = results%hessian
            end if

            !Calculates gradient in normal coordinates
            NMGradients = matmul(transpose(modesInRange), reshape(results%gradients,([3*self%chemicalSystem%nAtoms,1])))

            if (master()) then
               call kfHistory%Write('Coords', self%chemicalSystem%xyzAtoms)
               call kfHistory%Write('Energy', results%energy)
               call kfHistory%Write('Gradients', results%gradients)
               call kfHistory%Write('maxGrad', maxval(abs(results%gradients)))
               call kfHistory%Write('rmsGrad', sqrt( sum(results%gradients**2) / real(size(results%gradients),KREAL) ))
               call kfHistory%Write('nSelecNormalModes', size(modesInRange(1,:)))
               call kfHistory%Write('NMGradients', NMGradients)
               call kfHistory%FinishEntry()
            endif

            write (iuout,'(A,I4,A,f16.8,2(A,f12.8))') 'Step: ', self%iStep, ', Energy = ', results%energy, ', MaxNMGrad = ', &
                                                      maxval(abs(NMGradients)), ', MaxGrad = ', maxval(abs(results%gradients))

            !Check convergence and adjust the geometry if we are not converged yet.
            converged = all(abs(NMGradients) < self%gradConv)

            if (converged) then
               write (iuout, '(A)') 'Geometry converged!'
               write (iuout, '(A)')
               write (iuout,'(A,I5)') 'Number of normal modes (last step) found in the selected range', size(modesInRange(1,:))
            else
               !Generates a new step
               call newStep(self, cartHessian, NMGradients, modesInRange)
               deallocate(modesInRange) 
               deallocate(NMGradients)
            endif

            !Do we still have results from the previous iteration? If so, delete them now.
            !(At any time we want to have at most two results objects.)
            if (allocated(prevResults)) call prevResults%Delete()
            !Move the results from this iteration to prevResults, so that we have results free again for the next cycle.
            call move_alloc(results, prevResults)

         enddo

         !Clean up: Delete the last remaining results object
         call prevResults%Delete()
      endif

      deallocate(cartHessian)
      
      !Master closes the history and writes the optimized system to the ams.rkf file.
      if (master()) then
         call kfHistory%Close()
         call self%chemicalSystem%WriteToFile(self%kfUnit)
      endif

      call self%chemicalSystem%Print(iuout)
   
      call msg('Done!')

   end subroutine


   subroutine DeleteDerived(self)
      class(GeoOptNormalCoordsApplicationType), intent(inout) :: self

      ! Nothing to do here since we don't own any resources
      DontWarnUnused(self)

   end subroutine


   subroutine initHessian(self, cartHessian)

      !This code was copied from the SetUpInitialHessian subroutine in the OldGeoOptimizer.f90 file.

      use InternalCoordsModule

      class(GeoOptNormalCoordsApplicationType), intent(inout)  :: self
      integer(KINT)                                      :: numInternalCoords, iCoord
      real(KREAL)         , allocatable                  :: internalHessianDiag(:) 
      real(KREAL)         , allocatable , intent(out)    :: cartHessian(:,:)
      real(KREAL)                       , pointer        :: BMatrix(:,:)
      type(InternalCoordsType)                           :: internalCoordSet

      call New(internalCoordSet, self%chemicalSystem%xyzAtoms, self%chemicalSystem%GetAllAtomicNumbers(), self%chemicalSystem%GetAllMasses())

      numInternalCoords = GetNumOptimCoords(internalCoordSet)

      allocate(internalHessianDiag(numInternalCoords))
      allocate(BMatrix(numInternalCoords,3*self%chemicalSystem%nAtoms))

      call CalcApproxDiagHessian(internalCoordSet, internalHessianDiag)

      internalHessianDiag = sqrt(abs(internalHessianDiag))

      call CalcBMatrix(internalCoordSet, BMatrix)

      do iCoord = 1, 3*self%chemicalSystem%nAtoms
         BMatrix(1:numInternalCoords,iCoord) = BMatrix(1:numInternalCoords,iCoord) * internalHessianDiag(1:numInternalCoords)
      end do

      call Delete(internalCoordSet)
      deallocate(internalHessianDiag)

      cartHessian = matmul(transpose(BMatrix),BMatrix)

      deallocate(BMatrix)

   end subroutine


   subroutine modeSetup (self, modes, modesInRange)

      !Selects modes that fall within a frequency range interval selected by the user

      class(GeoOptNormalCoordsApplicationType), intent(inout)  :: self
      integer(KINT)                                      :: nSelecModes, iCounter, iFreq
      real(KREAL)         , allocatable , intent(out)    :: modesInRange(:,:)
      type(NormalModesType)             , intent(in)     :: modes
      
      !Number of modes within selected interval
      nSelecModes = count(modes%vibrFrequencies(:) > self%lowFreq .and. modes%vibrFrequencies(:) < self%highFreq)

      !modesInRange: S_{3n_{atoms} \times k_modes}
      allocate(modesInRange(3*self%chemicalSystem%nAtoms, nSelecModes))

      iCounter = 1
      do iFreq = 1, modes%nVibrModes
         if (modes%vibrFrequencies(iFreq) > self%lowFreq .and. modes%vibrFrequencies(iFreq) < self%highFreq) then
            modesInRange(:,iCounter) = modes%vibrModes(:,iFreq)
            iCounter = iCounter + 1
         endif
      enddo

   end subroutine


   subroutine newStep (self, cartHessian, NMGradients, modesInRange)

      !Produce a new step using the quadratic dependence and its RFO extension

      class(GeoOptNormalCoordsApplicationType), intent(inout)  :: self
      integer(KINT)                                      :: iCoord, nSelecModes
      real(KREAL)                       , intent(in)     :: CartHessian(:,:), NMGradients(:,:), modesInRange(:,:) 
      real(KREAL)         , allocatable                  :: RFO_denominator(:), QHessian(:,:), displCoords(:,:)
      real(KREAL)                                        :: Qii

      !Calculates the Q Hessian
      !QHessian: Q = S^t F S; [nModes x nModes] = [nModes x 3*nAtoms] * [3*nAtoms x 3*nAtoms] * [3*nAtoms * nModes]
      QHessian = matmul( matmul( transpose(modesInRange), cartHessian), modesInRange )

      !number of modes within selected range
      nSelecModes = size(QHessian(1,:))

      !Eqn. 7 (Step 5) in Referece [1]
      !displCoords: dq^{(i+1)}; [nModes x 1]
      allocate(displCoords(nSelecModes, 1)) 
      allocate(RFO_denominator(nSelecModes))

      do iCoord = 1, nSelecModes
         if (QHessian(iCoord, iCoord) < 1.0E-6_KREAL) call msg("WARNING: (very small) Qii element changed to 1.0E-6.")
         !Limits size of denominator so step sizes don't blow up
         Qii = max(abs(QHessian(iCoord,iCoord)), 1.0E-6_KREAL)
         !Q_{ii} + \sqrt{ Q_{ii}^2 + 4*g^2 }
         RFO_denominator(iCoord) = Qii + sqrt( (Qii)**2 + 4*(NMGradients(iCoord,1))**2 )
         !2*g / RFO_denominator
         !The negative sign in Reference [1] is not correct
         displCoords(iCoord,1) = 2 * (NMGradients(iCoord,1)) / RFO_denominator(iCoord) 
      enddo

      deallocate(QHessian)
      deallocate(RFO_denominator)

      !x^{(i+1)} = x^{(i)} - S dq{(i+1)}; [3 x nAtoms] = [3 x nAtoms] - [3*nAtoms x nModes] * [nModes x 1]
      !Eqn. 8 (Step 6) in Reference [1]
      self%chemicalSystem%xyzAtoms = self%chemicalSystem%xyzAtoms - reshape(matmul(modesInRange,displCoords), ([3,self%chemicalSystem%nAtoms]))

      deallocate(displCoords)

   end subroutine


   subroutine BFGSupdate(self, prevResults, currResults, cartHessian)

      !Broyden–Fletcher–Goldfarb–Shanno (BFGS) Hessian update algorithm is applied to generate new modes
      !cartHessian: f(i+1); [3nAtoms x 3nAtoms]
      !Eqn. 5 (Step 3) in Reference [1] is wrong, thus we got it from another reference.

      class(GeoOptNormalCoordsApplicationType) , intent(inout)  :: self
      class(AMSResultsType)            , intent(in)     :: prevResults, currResults
      real(KREAL)                      , intent(inout)  :: cartHessian(:,:)
      real(KREAL)         , allocatable                 :: deltaCoords(:,:), deltaGradients(:,:)
      real(KREAL)                                       :: BFGSdenom1(1,1), BFGSdenom2(1,1)

      !deltaCoords: dx^{(i)} = x^{(i)} - x^{(i-1)}; [3*nAtoms x 1] = [nAtoms, 3] - [nAtoms, 3]
      deltaCoords    = reshape((self%chemicalSystem%xyzAtoms - prevResults%chemicalSystem%xyzAtoms), ([3*self%chemicalSystem%nAtoms,1]))
      !deltaGradients: \Delta g^{(i)} = g^{(i)}_c - g^{(i-1)}_c; [3*nAtoms x 1] = [nAtoms, 3] - [nAtoms, 3]
      deltaGradients = reshape((currResults%gradients - prevResults%gradients), ([3*self%chemicalSystem%nAtoms,1]))

      !First we calculate the denominators in the BFGS equation
      !\Delta g^t dx; [1 x 1] =  [3*nAtoms x 1]^t * [3*nAtoms x 1]
      BFGSdenom1 = (matmul(transpose(deltaGradients), deltaCoords))
      !dx^t f dx; [1 x 1] = [3*nAtoms x 1]^t * [3*nAtoms x 3*nAtoms] * [3*nAtoms x 1]
      BFGSdenom2 = (matmul( matmul( transpose(deltaCoords), cartHessian ), deltaCoords))

      !If one of the denominators is too small, the Hessian is not updated.
      if (BFGSdenom1(1,1) > 1e-10 .or. BFGSdenom2(1,1) > 1e-10) then
         !The second term of the BFGS equation is then summed to the Cartesian Hessian from the previous step.
         !f = f - (f dx dx^t f) / BFGS_denominator
         ![3*nAtoms x 3*nAtoms] = [3*nAtoms x 3*nAtoms] - ([3*nAtoms x 3*nAtoms] * [3*nAtoms x 1] * [3*nAtoms x 1]^t * [3*nAtoms x 3*nAtoms]) / [1 x 1]
         cartHessian = cartHessian - (( matmul( matmul( matmul( cartHessian, deltaCoords ), transpose(deltaCoords) ), cartHessian ) ) / BFGSdenom2(1,1))
         !The first term of the BFGS equation is summed to the cartHessian (this order is done on purpose as the second term depends on the Hessian)
         !f = f + (\Delta g \Delta g^t) / BFGS_denominator; [3*nAtoms x 3*nAtoms] = [3*nAtoms x 3*nAtoms] + [3*nAtoms x 1] * [3*nAtoms x 1]^t / [1 x 1]
         cartHessian = cartHessian + ((matmul(deltaGradients, transpose(deltaGradients))) / BFGSdenom1(1,1))
      else
         call msg("WARNING: BFGS not updated in Step "//ftlString(self%iStep)//" as denominator was too small.")
      endif

      deallocate(deltaCoords, deltaGradients)

   end subroutine


   subroutine loadHessian(self, cartHessian)

      class(GeoOptNormalCoordsApplicationType), intent(inout)   :: self
      integer(KINT)                                       :: iuKF             ! File unit for filename RKF file
      logical                           , external        :: master
      real(KREAL)        , allocatable  , intent(out)     :: cartHessian(:,:)


      allocate(cartHessian(3*self%chemicalSystem%nAtoms,3*self%chemicalSystem%nAtoms))
      cartHessian = 0.0_KREAL

      if (.not. kfexfl(self%fileName%raw)) then   ! check if kf exists
         call StopIt('File '//self%fileName//' does not exist.')
      endif
      if (.not. isKFFile(self%fileName%raw)) then ! check if it is a KF file
         call StopIt('File '//self%fileName//' is not a KF file.')
      endif

      call kfopfl (iuKF, self%fileName%raw)
      if (.not.kfexvr(iuKF, 'AMSResults%Hessian')) call StopIt('File '//self%fileName//' does not contain a Hessian.')
      call kfread (iuKF, 'AMSResults%Hessian', cartHessian)
      call kfclfl (iuKF)

   end subroutine


end module
