character(100) function sweep_blanks(in_str)
character(100), intent(in) :: in_str
character(100) :: out_str
character :: ch
integer :: j

out_str = " "
    do j=1, len_trim(in_str)
!      get j-th char
      ch = in_str(j:j)
      if (ch .ne. " ") then
        out_str = trim(out_str) // ch
      endif
      sweep_blanks = out_str
    end do
  end function sweep_blanks