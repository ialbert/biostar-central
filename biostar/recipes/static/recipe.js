/*
Code that handles recipe interface goes here.
*/

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

    json_edit.setSize(null, 900);

    return json_edit
}

$(document).ready(function () {

    var template_edit = prepare_script();
    var json_edit = prepare_interface();
    template_edit.refresh();
    json_edit.refresh();

});