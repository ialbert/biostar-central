
// comment add
function show_add_comment(parent, post_id){
    csrf_html = parent.find("input[name='csrfmiddlewaretoken']").parent().html()
    $(".comment-form").remove();
    parent.append('\
    <form action="/new/comment/' + post_id + '/" method="post" class="comment-form" id="comment-form">            \
        ' + csrf_html + ' \
        <div class="controls"><textarea class="input-xlarge span8" id="textarea" name="content" rows="3"></textarea></div> \
        <div><a class="btn btn-success" href=\'javascript:document.forms["comment-form"].submit()\'><i class="icon-comment"></i> Add comment</a>          \
        (<span> markdown ok</span>) <a class="btn btn-warning" onclick="javascript:obj=$(\'.comment-form\').remove();"><i class="icon-remove"></i> Cancel</a>   </div>       \
    </form>            \
    '
    )
}

//comment deletion
function comment_delete(link){
    pid  = link.attr('target') // post id encoded in link
    par  = $('#'+pid)         // parent div
    body = par.children('div[name="content"]')
    $.post('/comment/delete/' + pid +'/',
        function(data){
            popover(link.parent(), data.msg, data.status);
            if (data.msg == 'Post destroyed') {
                par.hide('fast')
            } else {
                body.html('[ content deleted ]')
            }
        }, 'json');
}

function remove(elem, value){
    if (elem.value == value) {
        elem.value = ''
    }
    console.log(elem)
}
