# Understanding Recipes

## What is a recipe

A recipe consists of a "JSON definition file" and a "script template".

The simplest JSON definition file is

    {}

A simple script template might contain just:

    echo 'Hello World!'

## Recipe execution

Before executing the recipe the script template is rendered with the JSON data and is filled into the template.

    template + JSON -> script

The script is then executed at the command line.

## Interface file definition

The JSON definition file lists the parameters and allows the interface to be rendered.
Here is an example JSON definition file:

```
{
  foo: {
    label: Enter the name
    help: The name to appear after the greeting
    display: TEXTBOX
    value: World!
  }
}
```

the parameter name is `foo`, the default value is `World!`. The `display` field specifies the type of the HTML widget, the `label` and  `help` fields describe the interface. The interface generated from this specification file looks like this:

![Generated interface](interface1.png)

## Recipe template

A recipe is a script that has template markers for filling in parameters. In the case for the `foo` variable above, we can access its value via:

    echo 'Hello {{foo.value}}'

Recipes are using [Django templates][templates] and may contain Django template specific constructs.

## Recipe runtime

When the recipe is run the template will be substituted according to the interface value entered by the user. If the default value is kept it will produce the script:

    echo 'Hello World!'

[template]: https://docs.djangoproject.com/en/2.2/topics/templates/

## Results directory

Once the recipe runs a results directory is created that contains the following:

- the code for the recipe
- the standard out and error stream content
- all files created by the recipe

The results directory is a snapshot of all files generated when the recipe has been run, including the recipe itself.

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
