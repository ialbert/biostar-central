# `recipe` concepts

## What is a recipe

A recipe consists of a JSON definition file and a script template.

Both files may be empty.

## Recipe execution

Before executing the recipe the script template is rendered with the JSON data.

    template + JSON -> script

The script is then executed at the command line.

## JSON definition

The JSON definition file lists the paramters and allows the interface to be rendered. Here is an example JSON definition file:

    {
        name: {
        }
    }


## Recipe template

A recipe is a script that has template markers for filling in parameters. For example:

    echo "Hello {{name}}"

Recipes are using [Django templates][templates] and may contain Django template specific constructs.

[template]: https://docs.djangoproject.com/en/2.2/topics/templates/

## Data representation

A "data" unit in the `recipes` app is a directory that may contain one or more (any number of files).

## Data value

Each recipe parameter will have an automatic attribute called `value` that contains either the selected value (if  the parameter is user supplied) or the first file from the `table-of-contents`. For data consisting of a single file one may use the value directly.

    fastqc {{reads.value}}

## Data table-of-contents

Each recipe parameter will have an automatically generated attribute called `toc` (table of contents) that returns the list of the file paths in the data.

The file paths are absolute paths. The `toc` can be used to automate the processing of data. For example
a data directory named `reads` contains several FASTQ files with `.fq` extensions. To run `fastqc` on each file that matches that
the recipe may use:

    cat {{reads.toc}} | grep .fq | parallel fastqc {}

## Data source

When a recipe parameter indicates the source of the parameter as `PROJECT` it will be populated from the data in the project that matches the type.

    reference: {
        label: Reference Genome
        display: DROPDOWN
        type: FASTA
        source: PROJECT
    }

Only data that matches the tage `FASTA` will be shown in the dropdown menu.

## Data types

Data types are labels (tags) attached to each data that help filtering them in dropdown menus. More than one data type may be listed as comma separated values.
The data types may be any word (though using well recognized names: BED, GFF is recommended).

## File storage

Data that exists on a filesystem may be linked into the Biostar Engine from the command line. This means that no copying/moving of data is required. The only limitation is that of the filesystem.

## User permissions

Users may have read and write access to projects. A write access means that users may modify information, upload data and execute recipes.
staff and admin users can edit the recipe code.
