from django import forms
from thapp.models import Theory

class TheoryForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(forms.ModelForm, self).__init__(*args, **kwargs)
        for field_name, field in self.fields.items():
            field.widget.attrs['class'] = 'form-control'
    
    class Meta:
        model = Theory
        fields = '__all__'
        exclude = ['id']
