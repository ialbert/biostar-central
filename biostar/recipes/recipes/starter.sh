# Demonstration recipe

# Assign interface parameter value to bash variable.
COUNT={{count.value}}

# From here on down this is a regular bash script.

# Stop script on errors.
set -ue

# Repeat the message.
for n in `seq 10`; do
  echo "$n: Hello World!";
done