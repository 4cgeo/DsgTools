<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_estrut_transporte"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_org_comerc_serv"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="operacional">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="situacaofisica">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Planejada" value="1"/>
        <value key="Constru�da" value="2"/>
        <value key="Abandonada" value="3"/>
        <value key="Destru�da" value="4"/>
        <value key="Em constru��o" value="5"/>
        <value key="Constru�da, mas em obras" value="6"/>
        <value key="Desconhecida" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="matconstr">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Alvenaria" value="2"/>
        <value key="Madeira" value="6"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="turistica">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="cultura">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoedifcomercserv">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Farm�cia" value="1"/>
        <value key="Oficina mec�nica" value="2"/>
        <value key="Loja de materiais de constru��o e/ou ferragem" value="3"/>
        <value key="Centro comercial" value="4"/>
        <value key="Loja de conveni�ncia" value="5"/>
        <value key="Centro de conven��es" value="6"/>
        <value key="Motel" value="7"/>
        <value key="Loja de m�veis" value="8"/>
        <value key="Supermercado" value="9"/>
        <value key="Centro de exposi��es" value="10"/>
        <value key="Posto de combust�vel" value="11"/>
        <value key="Loja de roupas e/ou tecidos" value="12"/>
        <value key="Mercado p�blico" value="13"/>
        <value key="Quiosque" value="14"/>
        <value key="Quitanda" value="15"/>
        <value key="Com�rcio de carnes" value="16"/>
        <value key="Hotel" value="17"/>
        <value key="Banca de jornal" value="18"/>
        <value key="Venda de ve�culos" value="19"/>
        <value key="Banco" value="20"/>
        <value key="Pousada" value="21"/>
        <value key="Outros com�rcios" value="22"/>
        <value key="Outros servi�os" value="23"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="finalidade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Outros" value="99"/>
        <value key="Uso restrito" value="2"/>
        <value key="Residencial" value="3"/>
        <value key="Comercial" value="4"/>
        <value key="Servi�o" value="5"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>