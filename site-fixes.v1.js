/* V-Able menu fix v1. Only applies to vable.erdogan.edu.tr. */
(function(){'use strict';
if(location.hostname!=='vable.erdogan.edu.tr'||/^\/CP(?:\/|$)/i.test(location.pathname))return;
function install(){
var menu=document.getElementById('dikey-menu');if(!menu)return;
if(menu.classList.contains('vable-menu-fixed')){sync();return;}
var groups=Array.from(menu.children).filter(function(li){return li.querySelector(':scope > a.has-sub')&&li.querySelector(':scope > ul.sub-menu');});
var initiallyOpen=groups.find(function(li){return li.querySelector(':scope > a.has-sub').getAttribute('aria-expanded')==='true';});
if(!document.getElementById('vable-menu-fix-style')){var style=document.createElement('style');style.id='vable-menu-fix-style';style.textContent='#dikey-menu.vable-menu-fixed > li > ul.sub-menu{display:none!important}#dikey-menu.vable-menu-fixed > li.vable-menu-open > ul.sub-menu{display:block!important}';document.head.appendChild(style);}
menu.classList.add('vable-menu-fixed');
groups.forEach(function(li){
var sub=li.querySelector(':scope > ul.sub-menu');var duplicate=sub.querySelector(':scope > li.menu-ust-hedef');
if(duplicate){var generated=duplicate.querySelector('a');if(generated&&Array.from(sub.children).some(function(other){var a=other.querySelector('a');return other!==duplicate&&a&&a.href===generated.href;}))duplicate.remove();}
});
function setOpen(target){groups.forEach(function(li){li.classList.toggle('vable-menu-open',li===target);});sync();}
function sync(){Array.from(menu.children).forEach(function(li){var a=li.querySelector(':scope > a.has-sub');var sub=li.querySelector(':scope > ul.sub-menu');if(!a||!sub)return;var open=li.classList.contains('vable-menu-open');a.setAttribute('aria-expanded',String(open));sub.hidden=!open;});}
menu.addEventListener('click',function(event){var a=event.target.closest('a.has-sub');if(!a||a.parentElement.parentElement!==menu)return;event.preventDefault();event.stopImmediatePropagation();setOpen(a.parentElement.classList.contains('vable-menu-open')?null:a.parentElement);},true);
menu.addEventListener('keydown',function(event){var a=event.target.closest('a.has-sub');if(event.key==='Escape'){var open=groups.find(function(li){return li.classList.contains('vable-menu-open');});if(open){event.preventDefault();event.stopImmediatePropagation();setOpen(null);open.querySelector(':scope > a.has-sub').focus();}return;}if(a&&a.parentElement.parentElement===menu&&(event.key==='Enter'||event.key===' ')){event.preventDefault();event.stopImmediatePropagation();setOpen(a.parentElement.classList.contains('vable-menu-open')?null:a.parentElement);}},true);
setOpen(initiallyOpen||null);
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',install,{once:true});else install();
window.addEventListener('load',function(){requestAnimationFrame(install);},{once:true});
})();
