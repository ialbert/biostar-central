#
# A starter recipe with examples.
#

#
# You can fill in shell variables
#
READLEN={{readlen.value}}

echo "Read length: $READLEN"

#
# Substitute into content
#
echo "Referene genome: {{reference.value}}"

#
# But you may also use Django Template constructs.
#
{% if instrument.value == 'pacbio' %}

    echo "Yes, it is Pacific Biosciences!"

{% else %}

    echo "No, it is not Pacific Biosciences!"

{% endif %}

#
# Generate a table of content with all files in the job directory.
#
find . -name '*' > files.txt

#
# Print the contents to the screen
#
echo "****** File List: files.txt ****"
cat files.txt

# Make a nested directory
mkdir -p foo/bar
find . -name '*' > foo/bar/all.txt