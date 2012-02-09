function show_add_comment(parent, post_id){
    csrf_html = parent.find("input[name='csrfmiddlewaretoken']").parent().html()
    parent.html('\
    <form action="/new/comment/' + post_id + '/" method="post" id="comment-form">            \
        ' + csrf_html + ' \
        <div class="controls"><textarea class="input-xlarge span8" id="textarea" name="content" rows="3"></textarea></div> \
        <div><a class="btn btn-success" href=\'javascript:document.forms["comment-form"].submit()\'><i class="icon-comment"></i> Add comment</a>          \
        (<span> markdown ok</span>) <a class="btn btn-warning" href=\'javascript:$("#comment-form")\'><i class="icon-remove"></i> Cancel</a>   </div>       \
    </form>            \
    '    
    )
}