# some functions for dealing with strings
# there ought to be a regular expression way to doing this
function findWord(str::AbstractString, substr::AbstractString)
# returns the entire word containing a substring

substr_len = length(substr)
end_index = 0  # location of end of substring
start_index = 0  # starting index of word
for i=1:(length(str) - substr_len + 1)
  tmp = str[i:(i + substr_len - 1)]  # get the current substring

  if tmp == substr
    end_index = i + substr_len - 1
    start_index = getWordStart(str, end_index, ' ') + 1
    end_index = getWordStart(str, end_index, '/') - 1  # remove libmds.so from end
  end
end

return str[start_index:end_index]

end


function getWordStart(str::AbstractString, end_index::Integer, delim::Char)
# look backwards in the string for a space

  start_index = 1  # return 1 if no space found
  for i=end_index:(-1):1
    char_i = str[i]
#    println("index = ", i, " , ", "char = ", char_i, " , val = ", Int32(char_i))
    if char_i == delim
      start_index = i
      break
    end

  end

  return start_index
end



