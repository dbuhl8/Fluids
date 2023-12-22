subroutine deallocate_uTCB(u,Temp,Chem,B)
  implicit none
  type(velocity) :: u
  type(field)    :: B
  type(buoyancy) :: Temp, Chem

  deallocate(u%phys,u%spec,u%curl)
#ifdef AB_BDF3
  deallocate(u%rhs) 
#endif


#ifdef TEMPERATURE_FIELD
  deallocate(Temp%phys,Temp%spec)
#ifdef AB_BDF3
  deallocate(Temp%rhs) 
#endif
#endif

#ifdef CHEMICAL_FIELD
  deallocate(Chem%phys,Chem%spec)
#ifdef AB_BDF3
  deallocate(Chem%rhs) 
#endif
#endif

#ifdef MAGNETIC
  deallocate(B%phys,B%spec,B%curl) 
#ifdef AB_BDF3
  deallocate(B%rhs)
#endif
#endif


end subroutine deallocate_uTCB
