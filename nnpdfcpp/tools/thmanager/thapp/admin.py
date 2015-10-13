from django.contrib import admin

# Register your models here.

from .models import Theory

class TheoryAdmin(admin.ModelAdmin):
    using = 'theory'
    list_display = ('id','comments')
    def get_queryset(self, request):
        return super(TheoryAdmin, self).get_queryset(request).using(self.using)


admin.site.register(Theory, TheoryAdmin)
