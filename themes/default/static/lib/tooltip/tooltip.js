(function($) {

  var tooltip;
  var arrow;
  var arrowWidth;
  var arrowHeight;
  var content;
  var win;

  function getState(el, options) {
    var s = {};
    var elementHeight = el.outerHeight();
    var elementWidth  = el.outerWidth();
    var offset = el.offset();
    s.height = tooltip.outerHeight(true);
    s.width  = tooltip.outerWidth(true);
    s.offset = {};
    s.offset.top = offset.top;
    s.offset.left = offset.left;
    s.offset.right = s.offset.left + elementWidth;
    s.offset.bottom = s.offset.top + elementHeight;
    s.offset.hCenter = s.offset.left + Math.floor(elementWidth / 2);
    s.offset.vCenter = s.offset.top + Math.floor(elementHeight / 2);
    s.css = {};
    s.on  = {};
    s.off = {};
    s.arrow = {};
    return s;
  }

  function checkBounds(s, direction, margin, slide) {
    var bound, alternate;
    margin = parseInt(margin);
    slide  = parseInt(slide);
    switch(direction) {
      case 'top':
        bound = win.scrollTop();
        if(s.offset.top - s.height - margin - slide < bound) alternate = 'bottom';
        s.on.top  = s.offset.top - s.height - margin;
        s.off.top = s.on.top + slide;
        s.css.top = s.on.top - slide;
        s.css.left = getCenter(s, true);
        break;
      case 'left':
        bound = win.scrollLeft();
        if(s.offset.left - s.width - margin - slide < bound) alternate = 'right';
        s.on.left  = s.offset.left - s.width - margin;
        s.off.left = s.on.left + slide;
        s.css.top  = getCenter(s, false);
        s.css.left = s.on.left - slide;
        break;
      case 'bottom':
        bound = win.scrollTop() + win.height();
        if(s.offset.bottom + s.height + margin + slide > bound) alternate = 'top';
        s.on.top   = s.offset.bottom + margin;
        s.off.top  = s.offset.bottom - slide + margin;
        s.css.top  = s.on.top + slide;
        s.css.left = getCenter(s, true);
        break;
      case 'right':
        bound = win.scrollLeft() + win.width();
        if(s.offset.right + s.width + margin + slide > bound) alternate = 'left';
        s.on.left  = s.offset.right + margin;
        s.off.left = s.on.left - slide;
        s.css.left = s.on.left + slide;
        s.css.top = getCenter(s, false);
        break;
    }
    if(alternate && !s.over) {
      s.over = true;
      checkBounds(s, alternate, margin, slide);
    } else {
      s.direction = direction;
      getArrowOffset(s, direction);
      checkSlide(s, direction);
    }
  }

  function checkSlide(s, dir) {
    var offset;
    if(dir == 'top' || dir == 'bottom') {
      offset = win.scrollLeft() - s.css.left + 5;
      if(offset > 0) {
        s.css.left += Math.abs(offset);
        s.arrow.left -= offset;
      }
      offset = (s.css.left + s.width) - (win.scrollLeft() + win.width()) + 5;
      if(offset > 0) {
        s.css.left -= Math.abs(offset);
        s.arrow.left += offset;
      }
    } else if(dir == 'left' || dir == 'right') {
      offset = win.scrollTop() - s.css.top + 5;
      if(offset > 0) {
        s.css.top += Math.abs(offset);
        s.arrow.top -= offset;
      }
      offset = (s.css.top + s.height) - (win.scrollTop() + win.height()) + 5;
      if(offset > 0) {
        s.css.top -= Math.abs(offset);
        s.arrow.top += offset;
      }
    }
  }

  function getArrowOffset(s, dir) {
    if(dir == 'left' || dir == 'right') {
      s.arrow.top = Math.floor((s.height / 2) - (arrowHeight / 2));
    } else {
      s.arrow.left = Math.floor((s.width / 2) - (arrowWidth / 2));
    }
    s.arrow[getInverseDirection(dir)] = -arrowHeight;
  }

  function getInverseDirection(dir) {
    switch(dir) {
      case 'top':    return 'bottom';
      case 'bottom': return 'top';
      case 'left':   return 'right';
      case 'right':  return 'left';
    }
  }

  function getCenter(s, horizontal) {
    if(horizontal) {
      return s.offset.hCenter + (-s.width / 2);
    } else {
      return s.offset.vCenter + (-s.height / 2);
    }
  }

  function animateTooltip(s, options, el, fn) {
    var color = getDefault('color', options, el, 'black');
    var duration = getDefault('duration', options, el, 150);
    tooltip.attr('class', color + ' ' + s.direction);
    tooltip.stop(true, true).css(s.css);
    arrow.attr('style', '').css(s.arrow);
    tooltip.animate(s.on, {
      duration: duration,
      queue: false,
      complete: fn
    });
    tooltip.fadeIn(duration);
  }

  function animateTooltipOut(s, options, el, fn) {
    var duration = getDefault('duration', options, el, 100);
    tooltip.animate(s.off, {
      duration: duration,
      queue: false,
      complete: fn
    });
    tooltip.fadeOut(duration);
  }

  function unescapeHTML(html) {
    if(/&/.test(html)) {
      html = $('<p/>').html(html).text();
    }
    return html;
  }

  function setContent(el, title) {
    var html;
    try {
      var ref = $(document.body).find(title);
    } catch(e) {
      // May throw a malfolmed selector error
    }
    if(ref && ref.length > 0) {
      html = ref.html();
    } else {
      html = unescapeHTML(title);
    }
    content.html(html);
  }

  function getDefault(name, options, el, defaultValue) {
    return or(options[name], el.data('tooltip-'+name), defaultValue);
  }

  function or() {
    for(var i = 0; i < arguments.length; i++) {
      if(arguments[i] !== undefined) {
        return arguments[i];
      }
    }
  }

  jQuery.fn.tooltip = function(options) {
    options = options || {};
    this.each(function() {
      var el = $(this);
      var title = el.attr('title');
      if(!title) return;
      var animating = false;
      var state;
      var timer;
      el.unbind('mouseenter').mouseenter(function() {
        var delay = getDefault('delay', options, el, 300);
        clearTimeout(timer);
        timer = setTimeout(function() {
          var margin    = getDefault('margin', options, el, 20);
          var slide     = getDefault('slide', options, el, 10);
          var direction = getDefault('direction', options, el, 'right');
          var t         = el.attr('title');
          if(t) {
            title = t;
          }
          el.removeAttr('title');
          setContent(el, options.html || title);
          state = getState(el, options);
          checkBounds(state, direction, margin, slide);
          animateTooltip(state, options, el, function() {
            animating = false;
          });
          animating = true;
        }, delay);
      });
      el.unbind('mouseleave').mouseleave(function() {
        clearTimeout(timer);
        if(!state) return;
        if(animating) {
          tooltip.fadeOut(100, function() {
            animating = false;
          });
        } else {
          animateTooltipOut(state, options, el, function() {
            animating = false;
          });
        }
        state = null;
        animating = true;
      });
    });
  };

  $(document).ready(function() {
    tooltip = $('<div id="tooltip" />').appendTo(document.body).css('position', 'absolute').hide();
    arrow   = $('<div class="arrow" />').appendTo(tooltip);
    content = $('<div class="content" />').appendTo(tooltip);
    win     = $(window);
    arrowWidth = arrow.width();
    arrowHeight = arrow.height();
    $('.tooltip').tooltip();
  });

})(jQuery);
