/*
Code that handles recipe interface goes here.
*/


function show_script() {
    $(".edit-recipe #script").show();

    $(".edit-recipe #metadata").hide();
    $(".edit-recipe #description").hide();
    $(".edit-recipe #interface").hide();


}


function show_interface(){
    $(".edit-recipe #interface").show();

    $(".edit-recipe #metadata").hide();
    $(".edit-recipe #description").hide();
    $(".edit-recipe #script").hide();


}

function show_metadata(){
    $(".edit-recipe #metadata").show();

    $(".edit-recipe #description").hide();
    $(".edit-recipe #script").hide();
    $(".edit-recipe #interface").hide();


}

function show_description(){
    $(".edit-recipe #description").show();

    $(".edit-recipe #script").hide();
    $(".edit-recipe #interface").hide();
    $(".edit-recipe #metadata").hide();
}

function prepare_script() {

    var template_edit = CodeMirror.fromTextArea(
        document.getElementById("template"),
        {
            lineNumbers: true,
            mode: 'shell',
        }
    );

    function updateTemplate() {
        template_edit.save();
    }

    template_edit.on('change', updateTemplate);

    template_edit.setSize(null, 900);

    return template_edit
}


function prepare_interface() {
    var json_edit = CodeMirror.fromTextArea(
        document.getElementById("json"), {
            lineNumbers: true,
            mode: "engine",
            autoRefresh: true,
        }
    );

    function updateJson() {
        json_edit.save();
    }

    json_edit.on('change', updateJson);
    //json_edit.focus();

    json_edit.setSize(null, 900);
    json_edit.refresh();

    return json_edit
}

$(document).ready(function () {


});