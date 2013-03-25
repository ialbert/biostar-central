$(document).ready(function(){
    $('.add-comment').each(function(){
        elem = $(this)

        elem.click(function(){
            show_add_comment($(this).parent(), $(this).parent().children('input').val()); 
        });
    });

});


function show_add_comment(parent, post_id){
    csrf_html = parent.find("input[name='csrfmiddlewaretoken']").parent().html()
    parent.html('\
    <form action="/new/comment/' + post_id + '/" method="post">            \
        ' + csrf_html + ' \
        <div><textarea name="content" rows="4" cols="80"></textarea></div>            \
        <div><input type="submit" value="Add Comment"> (<span> markdown ok</span>) </div>        \
    </form>            \
     \
    '    
    )
}
