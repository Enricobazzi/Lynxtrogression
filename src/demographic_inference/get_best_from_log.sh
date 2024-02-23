# input log
log=${1}
# basename of input log
log_basename=$(basename "${log}" .log)
# path of input log
log_dir=$(dirname "${log}")

# split log file into two files: one with the best log and one with the rest
# the file with the best log has to contain every line in the log file after the
# line containing the string "BEST"
# the file with the rest has to contain every line in the log file before the
# line containing the string "BEST"
# if the file does not contain the string "BEST", no file is created

# get line number of line containing the string "BEST"
line_num=$(grep -n "BEST" ${log} | cut -d ":" -f 1)

# split log file into two files
tail -n +${line_num} ${log} > ${log_dir}/${log_basename}.best