$(document).ready(function(){
    $('.vote').each(function(){
        elem = $(this)
        
        post_id = elem.children('input').val()

        
        up_button = elem.children('.vote-up')
        down_button = elem.children('.vote-down')
        up_button.click(function(){
            do_vote($(this), post_id, 0); 
        });
        down_button.click(function(){
            do_vote($(this), post_id, 1); 
        });
    });
});

function popover(parent, msg, cls){
    parent.append('<div></div>')
    elem = parent.children('div')
    elem.addClass('popover ' + cls)
    elem.text(msg)
    elem.delay(2000).fadeOut(500, function(){
        $(this).remove() 
    });
}

function do_vote(button, post, type){
    $.post('/vote/' , {post:post, type:type},
    function(data){
        popover(button, data.msg, data.status)  
    }, 'json');
}