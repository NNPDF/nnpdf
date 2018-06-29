from django.conf.urls import url
from thapp import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'add/$', views.addtheory, name='addtheory'),
]
