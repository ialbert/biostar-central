
function add_comment(elem) {

    // remove comment body if exists.
    $("#comment-row").remove();

    var post_uid = elem.attr('data-value');
    var project_uid = elem.attr('project-uid');
    var container = elem.parent().parent();
    var comment_url = elem.attr("comment-url");
    var csrf_html = jQuery("[name=csrfmiddlewaretoken]").val();

    container.after(`<div id="comment-row" class="ui basic segment inputcolor">
    <form id="comment-form" class="ui form" action=${comment_url}  method="post">
        <input type="hidden" name="parent_uid" id="parent_uid" value=${post_uid}/>
        <input type="hidden" name="project_uid" id="project_uid" value=${project_uid}/>
        <input type="hidden" name="csrfmiddlewaretoken" value=${csrf_html} />

        <div class="">
            <div id="wmd-button-bar-2"></div>
            <textarea class="wmd-input-2" id="wmd-input-2"  name="content" rows="6"></textarea>
        </div>
        <div>
            <button class="ui submit green button" type="submit">
                <i class="check icon"></i>Add Comment
            </button>
            <a class="ui orange right floated button" onclick="javascript:obj=$(\'#comment-row\').remove();">
            <i class="undo icon"></i> Cancel
            </a>
        </div>
    </form>
    </div>`
    );

    var converter = new Markdown.Converter();
    var editor = new Markdown.Editor(converter, '-2');
    editor.run();

};


$(document).ready(function () {


     $('select')
        .dropdown()
    ;

//    $(".items > .item").click(function (event) {
//        var obj = $(this).find("a:first");
//        if (typeof obj !== 'undefined') {
//            window.location = obj.attr("href");
//       }
//    });

    $(".copy-data").click(function (event) {

        event.preventDefault();
        var elem = $(this);
        var data_uid = elem.attr('data-uid');

        $.ajax("/data/copy", {
                type: 'GET',
                dataType: 'json',
                data: {data_uid: data_uid},
                success: function (data) {
                $("#copy-message-"+ data_uid).append(`
                <div class="ui basic segment">
                    <div class="ui fluid green message">
                    ${data.message}
                    </div>
                </div>

                `).fadeOut(5000);
                        },
                });

    });

    $(".add-comment").click(function (event) {
        add_comment($(this));
        });



});
