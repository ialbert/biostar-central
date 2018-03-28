## Concepts

### Distributed access

We envisions the Biostar Engine as a distributed system that runs on individual computers rather than one, single centralized location.

### Data Representation

A "data" in the Biostar Engine is a directory that may contain one or more (any number of files).

### Data Table of Contents

Each recipe parameter will have an automatic attribute called `toc` (table of contents) that returns the list of the file paths in the data.
The file paths are absolute paths. the `toc` can be used to automate the processing of data. For example
a data directory named `reads` contains FASTQ files with `.fq` extensions. To run `fastqc` on each file that matches that
the script may use:

    cat {{reads.toc}} | grep .fq | parallel fastqc {}

### Data Value

Each recipe parameter will have an automatic attribute called `value` that contains either the selected value (if  the parameter is user supplied) or
the first file from the `toc`. For data consisting of a single file one may use the value directly.

    fastqc {{reads.value}}

### Data Source

When a recipe parameter indicates the source of the parameter as `PROJECT` it will be populated from the data in the project that matches the type.

    reference: {
        label: Reference Genome
        display: DROPDOWN
        type: FASTA
        source: PROJECT
    }

Only data that matches the word `FASTA` will be shown in the dropdown menu.

### Data Types

Data types are labels (tags) attached to each data that help filtering them in dropdown menus. More than one data type may be listed as comma separated values.
The data types may be any word (though using well recognized names: BED, GFF is recommended).

### File Storage

Data that exists on a filesystem may be linked into the Biostar Engine from the command line. This means that no copying/moving of data is required.
The only limitation is that of the filesystem.

### User permissions

Users may have read and write access to projects. A write access means that users may modify information, upload data and execute recipes.

### Moderator permissions

Moderator permissions are site wide and specify the right to approve modifications to recipes. Because recipes execute arbitrary code only approved
recipes may be executed.

### Recipe diffs

Whenever a recipe changes it creates a "diff" (difference). The recipe will not be executable until a moderator reviews and accepts the diff.

### User types

Users may fall in different categories: VISITOR, RESEARCHER, MODERATOR etc.

Not all user types have access to all functionality. For example a VISITOR may see projects but not execute recipes.





