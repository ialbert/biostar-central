

function submit_comment (elem) {

    var comment_form = $("#comment-form");
    action_url = comment_form.attr("action");
    var container = $("#comment-row");

    $.ajax({
        url : action_url, // the endpoint
        type : "POST", // http method
        data : comment_form.serialize(), // data sent with the post request

        success : function(data) {
            // remove the value from the input
            alert(data.message);
            }


    });
        return false;
};



function add_comment(elem) {

    // remove comment body if exists.
    $("#comment-row").remove();

    var post_uid = elem.attr('data-value');
    var project_uid = elem.attr('project-uid');
    var container = $("#comment-container-"+ post_uid);
    var comment_url = elem.attr("comment-url");
    var redir_url = elem.attr("redir-url")
    var csrf_html = jQuery("[name=csrfmiddlewaretoken]").val();


    container.after(` <div id="comment-row" class="ui basic segment inputcolor">
    <form id="comment-form" class="ui form" action=${comment_url}  method="post">
        <input type="hidden" name="parent_uid" id="parent_uid" value=${post_uid} />
        <input type="hidden" name="project_uid" id="project_uid" value=${project_uid} />
        <input type="hidden" name="redir_url" id="redir_url" value=${redir_url} />
        <input type="hidden" name="csrfmiddlewaretoken" value=${csrf_html} />

        <div class="">
            <div id="wmd-button-bar-2"></div>
            <textarea class="wmd-input-2" id="wmd-input-2"  name="content" rows="6"></textarea>
        </div>
        <div>

            <a class="ui submit green button" onclick="return submit_comment($(this));">
                <i class="check icon"></i>Add Comment
            </a>
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
        var copy_url = elem.attr('copy-url');

        $.ajax(copy_url, {
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

                `).fadeOut(2000);
                        },
                });

    });

    $(".add-comment").click(function (event) {
        add_comment($(this));
        });

});
