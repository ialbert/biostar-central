$(document).ready(function(){
    $('.vote').each(function(){
        elem = $(this)
                
        up_button = elem.children('.vote-up')
        down_button = elem.children('.vote-down')
        up_button.click(function(){
            do_vote($(this), $(this).parent().children('input').val(), 0); 
        });
        down_button.click(function(){
            do_vote($(this), $(this).parent().children('input').val(), 1); 
        });
    });
});

function toggle_button(button){
    // Toggles the on/off status of a voting button
    if(button.hasClass('vote-on')){
        button.removeClass('vote-on');
        button.addClass('vote-off');
    }else if(button.hasClass('vote-off')){
        button.removeClass('vote-off');
        button.addClass('vote-on');
    }
}

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
    toggle_button(button) // Pre-emptitively toggle the button to provide feedback
    $.post('/vote/' , {post:post, type:type},
    function(data){
        popover(button.parent(), data.msg, data.status)
        if(data.status == 'error'){
            toggle_button(button) // Untoggle the button if there was an error
        } else {
            button.parent().children('.vote-count').text(data.new_score)
        }
    }, 'json');
}