
function toggle_css(elem, start, end){
    if(e.hasClass(start)){ e.removeClass(start); e.addClass(end); }
}

// modifies the votecount value
function mod_votecount(button, k){
    count = parseInt(button.siblings('.vote-count').text()) || 0
    count += k
    button.siblings('.vote-count').text(count)
}

function toggle_button(button){
    // Toggles the on/off status of a voting button
    if(button.hasClass('vote-on')){
        button.removeClass('vote-on');
        button.addClass('vote-off');
    } else if(button.hasClass('vote-off')){
        button.removeClass('vote-off');
        button.addClass('vote-on');
    }
    
    if (button.hasClass('vote-up') && button.hasClass('vote-on')){
        button.attr('data-original-title', 'Click to remove');
    }
    
    if (button.hasClass('vote-up') && button.hasClass('vote-off')){
        button.attr('data-original-title', 'Click to upvote');
    }
    
     if (button.hasClass('vote-bookmark') && button.hasClass('vote-on')){
        button.attr('data-original-title', 'Click to remove bookmark');
    }
    
    if (button.hasClass('vote-bookmark') && button.hasClass('vote-off')){
        button.attr('data-original-title', 'Click to add bookmark');
    }
    
    // Turn off opposite buttons if they're on
    if(button.hasClass('vote-on')){ 
        if(button.hasClass('vote-up')) 
            toggle_button(button.siblings('.vote-down.vote-on'))
        if(button.hasClass('vote-down'))
            toggle_button(button.siblings('.vote-up.vote-on'))       
    }
    // Update the vote counts immediately
    if(button.is('.vote-up.vote-on, .vote-down.vote-off'))
        mod_votecount(button, +1)
    if(button.is('.vote-up.vote-off, .vote-down.vote-on'))
        mod_votecount(button, -1)
    
}

function popover(parent, msg, cls){
    parent.append('<div></div>')
    elem = parent.children('div').last()
    elem.addClass('vote-popover ' + cls)
    elem.text(msg)
    elem.delay(1000).fadeOut(1000, function(){
        $(this).remove() 
    }); 
}

function do_vote(button, post, type){
    toggle_button(button) // Pre-emptitively toggle the button to provide feedback
    $.ajax('/vote/', {
        type: 'POST',
        dataType: 'json',
        data: {post:post, type:type},
        success: function(data){
            if(data.status == 'error'){ // Soft failure, like not logged in
                popover(button.parent(), data.msg, data.status) // Display popover only if there was an error
                toggle_button(button) // Untoggle the button if there was an error
            }
        },
        error: function(){ // Hard failure, like network error
            popover(button.parent(), 'Unable to submit vote', 'error');
            toggle_button(button);
        }
    });
}

