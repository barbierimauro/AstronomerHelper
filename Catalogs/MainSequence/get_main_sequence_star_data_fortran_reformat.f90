program main_sequence
implicit none
character(len=211) :: line
character(len=20) :: fields(32)
character(len=20) :: missing_data = '...'
integer :: i, j, unit_input, unit_output
character(len=*), parameter :: input_file = 'msdata.tmp'
character(len=*), parameter :: output_file = 'msdata.csv'
integer :: pos

open(newunit=unit_input, file=input_file, status='old', action='read')
open(newunit=unit_output, file=output_file, status='replace', action='write')

read(unit_input, '(A)') line
write(unit_output, '(A)') line

do
	read(unit_input, '(A)', iostat=i) line
	if (i /= 0) exit
	! Split the line into fields
	j = 1
	do i = 1, size(fields)
    	pos = scan(line, ' ')
    	if (pos == 0) pos = len_trim(line) + 1
    	fields(i) = adjustl(line(1:pos-1))
    	line = adjustl(line(pos:))
    	if (len_trim(line) == 0) exit
    	j = j + 1
	end do

	! Replace missing data with an empty field
	do i = 1, size(fields)
    	if (fields(i) == missing_data) fields(i) = ''
	end do

	! Write the processed fields to the CSV file (excluding the last column)
	write(unit_output, '(A)', advance='no') trim(fields(1))
	do i = 2, size(fields) - 1
    	write(unit_output, '(A)', advance='no') ',' // trim(fields(i))
	end do
	write(unit_output, '(A)') ''
end do

close(unit_input)
close(unit_output)
end program main_sequence
